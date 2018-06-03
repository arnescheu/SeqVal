__author__ = "Arne Scheu"
__date__ = "2018-06-04"
__version__ = "1.0"

import csv
import datetime
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib.lines import Line2D

# https://github.com/Cadair/jupyter_environment_kernels/issues/10
try:
    FileNotFoundError
except NameError:  # python2
    FileNotFoundError = IOError
    print("WARNING defined FileNotFoundError as IOError. Are you running python2? \n")

def master():
    for task in define_task(
            os.path.join(wd, "task.txt")):  # define_task: Defining which files to handle from task.txt file
        print("Processing task %s %s %s" % (task["gene"], task["clone"], task["plasmid"]))

        # Defining gene of interest based on genbank and region or sequence field.
        task["goi"], task["start"], task["record"] = define_goi(task)
        if len(task["goi"]) % 3 != 0:
            print("WARNING defined region not multiple of 3! Likely bad output, check %s" % str(task))

        # define the sequence traces belonging to gb file
        ab1_list = []
        for key in task:
            if "sequencing" in key:  # sequencing1, sequencing2,...
                if task[key] is not "":
                    split = task[key].split(";")  # sequencing filename1;filename2 also works
                    for path in split:
                        ab1_list.append(path)

        # defaulting to all ab1 in directory if none defined
        if not ab1_list:
            for files in os.walk(wd):
                if ".ab1" in files:
                    ab1_list.append(files)
            print("WARNING no .ab1 traces were defined for %s. DEFAULT to all traces - %s" % (
                task["genbank"], str(ab1_list)))

        # For each sequence trace, creating corresponding alignment in nblast
        # Then creates image of traces at mismatch positions
        alignment_list = []
        for ab1_path in ab1_list:
            print("\tProcessing trace %s" % ab1_path)
            if local_nblast(task["goi"], ab1_path) != 0:
                print("WARNING trace not found, skipping trace")
                continue
            alignment = interpret_blast()
            if not alignment:
                continue
            ab1 = SeqIO.read(ab1_path, 'abi')
            traces = handle_alignment(ab1, alignment, task["gene"])
            alignment_list.append({"ab1": ab1, "alignment": alignment, "traces": traces})

        # Create expasy, removing first residue ([3:]
        with open(os.path.join(wd, "seq_align", "temp", "nomet_expasy.html"), "w+") as ex:
            ex.write(expasy(task["goi"][3:].translate(to_stop=True)).decode("utf-8"))

        # Compiles all previous data
        task["user"] = username
        html_love(task, alignment_list)

    print("All tasks complete.")

def define_task(taskfile):
    task_list = []
    with open(taskfile, "r") as csvfile:
        task_reader = csv.DictReader(csvfile, delimiter="\t")
        # https://stackoverflow.com/questions/15865750/dictreader-change-keys-to-upper
        for row in task_reader:
            task_list.append(row)
    return task_list


def define_goi(task):  # not bothering with -1 strand
    if not task["genbank"] is "":
        try:
            record = SeqIO.read(task["genbank"], "genbank")
            for feature in record.features:
                if "GOI" in feature.qualifiers["label"]:  # define label name by genbank file?
                    start, end = int(feature.location.start), int(
                            feature.location.end)  # location.start/end are already adjusted to index, -1
                    return (record.seq[start:end], start, record)
            if not task["region"] is "":  # define sequence from region
                start, end = int(task["region"].split(":")[0]) - 1, int(task["region"].split(":")[1]) - 1
                return (record.seq[start:end], start, record)
            elif not task["sequence"] is None:  # define region from sequence
                goi = Seq(task["sequence"])
                start = int(str(record.seq.find(str(goi).upper())))
                if start == -1:
                    raise Exception("Sequence provided but not found in gb!")
                return (goi, start, record)
            else:
                print("WARNING gene not defined by region/sequence/annotation. Default to entire genbank region")
                return (record.seq, 0, record)
        except FileNotFoundError:
            if not task["sequence"] is None:
                print("WARNING %s not found! DEFAULT to sequence provided" % task["genbank"])
                return (Seq(task["sequence"]), 0, False)
            else:
                # TODO if permissive: return False
                raise FileNotFoundError(task["genbank"] + " not found!")
    elif not task["sequence"] is None:
        print("WARNING no genbank provided. DEFAULT to sequence provided")
        return (Seq(task["sequence"]), 0, False)
    else:
        raise Exception("Task lacking either genbank with region and sequence.")
        # return False


""" in construction. Might be easier to just use pairwise2
from Bio.Blast import NCBIWWW
def online_nblast(goi,ab1_path):
    ab1 = SeqIO.read(ab1_path, 'abi')
    FASTA = ">goi\n" + str(goi) + "\n>seq\n" + str(seq)
    print("Attempting online blast,",FASTA)
    result_handle = NCBIWWW.qblast(program="blastn","nt", FASTA)
    with open(os.path.join("seq_align","temp","blast_temp.txt"),"w") as fh:
        fh.write(result_handle.read())
    return True
"""


def local_nblast(goi, seqfile):  # better to switch to pairwise2 instead of nblast?
    query_temp = os.path.join("seq_align", "temp", "BLAST_query.fa")
    seq_temp = os.path.join("seq_align", "temp", "BLAST_sequence.fa")

    # converts goi sequence to file for blast
    with open(query_temp, "w+") as fh:
        fh.write(">query\n" + str(goi))

    # converts trace to file for blast
    try:
        SeqIO.convert(seqfile, "abi", seq_temp, "fasta")
    except FileNotFoundError:
        return seqfile + " not found. Skipping BLAST."

    # resultfile location
    out_temp = os.path.join("seq_align", "temp", "blast_temp.txt")
    # execute nblast using installation
    command = (
            'blastn -subject "%s" -query %s -dust no -outfmt 3 -num_alignments 1 >%s' % (
        seq_temp, query_temp, out_temp))  # output 7 creates xml biopython blast output. Biopython can read
    return os.system(command)


def interpret_blast(
        blastfile="blast_temp.txt"):  # converts blast output so subject/query. Very bad practice, could just use other modules
    with open(os.path.join(os.getcwd(), "seq_align", "temp", blastfile), "r") as raw_blast:
        start, space = False, 0
        subject = {"query": "", "sequence": "", "misalign": []}
        subject["tiny_blast"] = ""
        subject["seq_range"] = [0, 0]
        subject["que_range"] = [0, 0]
        for line in raw_blast:
            space += 1
            if "Query_1" in line:
                start, space = True, 0
                linesplit = list(filter(None, line.strip("\n").split(" ")))
                subject["query"] = subject["query"] + linesplit[2]
                subject["que_range"][1] = int(linesplit[3])
                if not subject["que_range"][0]:
                    subject["que_range"][0] = int(linesplit[1])
            if start:
                if space == 1:
                    linesplit = list(filter(None, line.strip("\n").split(" ")))
                    if not subject["seq_range"][0]:
                        subject["seq_range"][0] = int(linesplit[1])
                    subject["sequence"] = subject["sequence"] + linesplit[2]
                    subject["seq_range"][1] = int(linesplit[3])

                    # format blast to html coloring
                    linesplit_rec = ""
                    cont = False
                    for char in linesplit[2]:
                        if char != ".":
                            if not cont:
                                linesplit_rec = linesplit_rec + "#style#"
                                linesplit_rec = linesplit_rec + char
                                cont = True
                            elif cont:
                                linesplit_rec = linesplit_rec + char
                        else:
                            if cont:
                                linesplit_rec = linesplit_rec + "</font>"
                                linesplit_rec = linesplit_rec + char
                                cont = False
                            else:
                                linesplit_rec = linesplit_rec + char
                    if cont == True:
                        linesplit_rec = linesplit_rec + "</font>"  # Bugfix character at end of line didn't close color
                    subject["tiny_blast"] = subject["tiny_blast"] + line.replace(linesplit[2], linesplit_rec)

                elif space > 2:  # stop after last query:subject block
                    break
                else:
                    if "Subject" in line:
                        print("WARNING aligns multiple times in BLAST. Skipping trace")
                        return False
                    subject["tiny_blast"] = subject["tiny_blast"] + line

        if len(subject["query"]) != len(subject["sequence"]):
            print("WARNING discontinuous alignment. Does trace belong to genbank? Skipping trace")
            return False
        for i, character in enumerate(subject["sequence"]):
            if character != ".":
                subject["misalign"].append(i)
        style = "<font color='red'>"
        subject["html_tiny_blast"] = subject["tiny_blast"].replace("\n", "<br>").replace(" ", "&nbsp;").replace(
                "#style#",
                style)

    if subject["seq_range"][0] < subject["seq_range"][1]:
        subject["orientation"] = 1
    else:
        subject["orientation"] = -1

    return subject


def handle_alignment(ab1, alignment, taskname):
    traces = []
    guide = ""
    for i in range(0, len(alignment["sequence"])):
        if alignment["sequence"][i] == ".":
            guide = guide + "."
        elif alignment["sequence"][i] == "-":
            guide = guide + "-"
        elif alignment["query"][i] == "-":
            guide = guide + "+"
        elif alignment["sequence"][i] != ".":
            guide = guide + "."

    trace_offset = 0
    query_offset = 0

    for position, char in enumerate(guide):
        if position in alignment["misalign"]:
            traces.append(chromatogram(position, trace_offset, query_offset, ab1, alignment, taskname))
        if char == ".":
            pass
        elif char == "+":
            query_offset = query_offset - 1
        elif char == "-":
            trace_offset = trace_offset - 1

    return traces


def chromatogram(blast_position, offset, query_offset, ab1, alignment, taskname):
    relative_position = blast_position + offset
    start = alignment["seq_range"][0] - 1
    mismatch_number = alignment["que_range"][0] + blast_position + query_offset
    if alignment["orientation"] == 1:
        absolute_position = start + relative_position  # base start-1 = array start
    else:
        absolute_position = start - relative_position  # base start-1 = array start
    readposition = ab1.annotations["abif_raw"]["PLOC1"][absolute_position]
    mutation = alignment["query"][blast_position] + str(mismatch_number) + alignment["sequence"][blast_position]

    """
    print("actual position and read", str(absolute_position), ab1.seq[absolute_position - 2:absolute_position + 2],
          ab1.seq[absolute_position])  # 1227,1224
    print("blast_position and query, subject", blast_position, alignment["query"][blast_position],
          alignment["sequence"][blast_position])
    print("read position", readposition)
    """

    # https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    plt.rc('font', size=14)  # controls default text sizes
    plt.rc('axes', titlesize=14)  # fontsize of the axes title
    plt.rc('axes', labelsize=13)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    plt.rc('legend', fontsize=12)  # legend fontsize
    plt.rc('figure', titlesize=16)  # fontsize of the figure title

    # http://biopython.org/wiki/ABI_traces
    channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']
    trace = defaultdict(list)
    for c in channels:
        trace[c] = ab1.annotations['abif_raw'][c]
    clr = {"DATA9": "black", "DATA10": "green", "DATA11": "red", "DATA12": "blue"}  # G,A,T,C

    # TODO define view area a bit more flexible. Hardcoded numbers right now
    ymax = 0  # autoscaling
    for channel in "DATA9", "DATA10", "DATA11", "DATA12":  # Not sure why I did this instead of plotting directly
        ymax2 = max(ab1.annotations['abif_raw'][channel][readposition - 200:readposition + 200])
        if ymax2 > ymax:
            ymax = ymax2  # Maximum channel value
        plt.plot(ab1.annotations['abif_raw'][channel], color=clr[channel])
    ymax = ymax * 1.125

    # setting view to 10 residues before and after
    xmin = ab1.annotations["abif_raw"]["PLOC1"][absolute_position - 10]
    try:
        xmax = ab1.annotations["abif_raw"]["PLOC1"][absolute_position + 10]
    except IndexError:
        print("WARNING reached end of chromatogram, DEFAULT to maximum. UNTESTED, please show to Arne.")
        xmax = len(ab1.annotations["abif_raw"]["DATA9"])

    # plt.title("Misaligned " + alignment["query"][blast_position] + str(mismatch_number) + alignment["sequence"][
    #    blast_position] + ", ab1" + str(
    #        readposition)  + " offset " + str(offset),  y=1.2)

    plt.title(alignment["query"][blast_position] + str(mismatch_number) + alignment["sequence"][blast_position], y=1.15)

    if alignment["orientation"] == 1:
        plt.xlim(xmin, xmax)  # For some reason, setting this again makes the ticks place correctly
    else:
        plt.xlim(xmax, xmin)

    ax1 = plt.gca()
    if alignment["orientation"] == 1:
        ax1.set_xlabel('Sequencing')
    else:
        ax1.set_xlabel('Sequencing reverse complement')

    plt.ylim(0, ymax)

    # Adding legend for G/A/T/C
    # stolen from somewhere on stackoverflow
    custom_lines = [Line2D([0], [0], color=clr["DATA10"], lw=3),
                    Line2D([0], [0], color=clr["DATA11"], lw=3),
                    Line2D([0], [0], color=clr["DATA12"], lw=3),
                    Line2D([0], [0], color=clr["DATA9"], lw=3), ]

    if alignment["orientation"] == 1:
        ax1.legend(custom_lines, ["A", "T", "C", "G"])
    else:
        ax1.legend(custom_lines, ["T", "A", "G", "C"])

    plt.axvline(x=readposition)

    ax2 = ax1.twiny()
    ax2.set_xlabel('Query')

    # Adding pointers on 2nd axis
    # https://stackoverflow.com/questions/27728312/matplotlib-ticks-not-align-with-data-point
    visible_ticks = []
    for r in range(absolute_position - 10, absolute_position + 11, 1):
        visible_ticks.append(int(ab1.annotations["abif_raw"]["PLOC1"][r]))
    ax2.set_xticks(visible_ticks)
    ax1.set_xticks(visible_ticks)

    fitquery = []  # This handles mismatches where sequence has a deletion, e.g. GA-A requires two characters in aligned query
    hold = ""
    alignment["seq_reconstituted"] = alignment["sequence"]
    for i, char in enumerate(alignment["sequence"]):
        if char == "-":
            hold = hold + alignment["query"][i]
        else:
            if hold != "" and alignment["seq_reconstituted"][
                i] == ".":  # This passage puts corresponding read oppossite of deletion if equal, e.g. AC-. to AC-C, but not AC-T to AC-C
                alignment["seq_reconstituted"] = alignment["seq_reconstituted"][:i] + alignment["query"][i] + alignment[
                                                                                                                  "seq_reconstituted"][
                                                                                                              i + 1:]
            hold = hold + alignment["query"][i]
            fitquery.append(hold)
            hold = ""

    alignment["seq_reconstituted"] = alignment["seq_reconstituted"].replace("-", "")

    if alignment["orientation"] == -1:
        ticks = alignment["seq_reconstituted"][relative_position - 10:relative_position + 11]
        alignment["seq_reconstituted"] = alignment["seq_reconstituted"][::-1]
        ticks = ticks[::-1]
        # A piece of control
        if ticks != alignment["seq_reconstituted"][len(alignment["seq_reconstituted"]) - relative_position - 11:len(
                alignment["seq_reconstituted"]) - relative_position + 10]:
            raise Exception("~~CONGRATULATIONS~~, this should never happen! Please talk to Arne")
    else:
        ticks = alignment["seq_reconstituted"][relative_position - 10:relative_position + 11]
    # if len(ticks)<21: print(ticks, "Short trace. Extending by", 21-len(ticks))
    while len(ticks) < 21:
        ticks = ticks + str(ab1.annotations["abif_raw"]["PBAS1"][absolute_position + (len(ticks)) - 10])
    if verbose:
        print("\t\tPosition %s chromatogram %s ab1 %s task %s" % (mutation.ljust(6, " "), ticks, taskname, ab1.name))

    """
    ticks2 = ""
    if alignment["orientation"] == 1:
        for i in range(-10, 11):
            ticks2 = ticks2 + str(ab1.annotations["abif_raw"]["PBAS1"][absolute_position + i])
    else:
        for i in range(-10, 11):
            pos = str(ab1.annotations["abif_raw"]["PBAS1"][absolute_position + i])
            pos = str(Seq(pos).complement()) #TODO seq class doesn't correctly handle W, B, ... but unused...
            ticks2 = ticks2 + pos
    print(ticks2, "Full trace")
    """

    ax1.set_xticklabels(ticks)

    if alignment["orientation"] == 1:
        ticks = fitquery[relative_position - 10:relative_position + 11]
        ax2.set_xticklabels(ticks)
        plt.xlim(xmin, xmax)
    else:
        ticks = fitquery[relative_position - 10:relative_position + 11][::-1]
        ax2.set_xticklabels(ticks)
        plt.xlim(xmax, xmin)

    plt.gcf().set_size_inches(6, 5)
    plt.tight_layout()

    # saving file
    figname = os.path.join("seq_align", "results", "traces", taskname + "_" + mutation + "_" + ab1.name + ".png")
    plt.savefig(os.path.join(os.getcwd(), figname))

    plt.gcf().clear()
    # plt.show()
    return {"figurefile": figname, "mutation": mutation,
            "position": mismatch_number}  # Bugfix blast_position to mismatch_number


def expasy(sequence):  # without initial methionine
    url = "https://web.expasy.org/cgi-bin/protparam/protparam"
    payload = {'sequence': sequence}
    r = requests.get(url, params=payload)
    return r.content


def html_love(html_dict, alignment_list):
    html_dict["date_seqval"] = "%s-%s-%s" % (datetime.datetime.now().year, str(datetime.datetime.now().month).zfill(2),
                                             str(datetime.datetime.now().day).zfill(2))
    goi = html_dict["goi"]  # unneccessary renaming. Can't be bothered right now
    html_dict["sequence"] = goi
    html_dict["len_sequence"] = len(goi)

    # Truncate expasy for replacing into template
    html_dict["expasy"] = ""
    with open(os.path.join(os.getcwd(), "seq_align", "temp", "nomet_expasy.html"), "r") as exh:
        body = False
        for line in exh:
            if "User-provided sequence:</h2>" in line:
                body = True
            if body:
                html_dict["expasy"] = html_dict["expasy"] + line
            if "<!-- sib_body -->" in line:
                html_dict["expasy"] = html_dict["expasy"] + "</div>"
                body = False

    # Translate features to HTML
    with open(os.path.join(os.getcwd(), "seq_align", "templates", "seqtemplate.txt"), "r") as fh:
        seqtemplate = fh.read()

    record = html_dict["record"]
    if record:
        sorted_features = sorted(record.features, key=lambda x: x.location.start)
        featurecomprehension = []
        for i in range(0, len(goi)):
            featurecomprehension.append([])
        filtered_features = []
        i = 0
        for feature in sorted_features:
            if not (feature.type == "CDS" or feature.type == "Protein"):
                continue
            if feature.location.strand == -1 or "GOI" in feature.qualifiers["label"]:
                continue
            featurestart = int(feature.location.start) - html_dict["start"]
            featureend = int(feature.location.end) - html_dict["start"]
            if featurestart > len(goi) or featureend > len(goi):
                continue
            for element in featurecomprehension[featurestart:featureend]:
                element.append(i)
            filtered_features.append(feature)
            i = i + 1

        # labeling of gene sequence
        sequence_buffer = ""
        annotation_buffer = ""
        level = [-1, -1]
        layer = 0
        col = ["", ""]
        for i, position in enumerate(featurecomprehension):
            if not labelcolor:
                if len(position) > 0:
                    col[0] = colors[position[0]]
                    if len(position) > 1:
                        col[1] = colors[position[1]]
            else:
                if len(position) > 0:
                    col[0] = filtered_features[position[0]].qualifiers["ApEinfo_revcolor"][0]
                    if len(position) > 1:
                        col[1] = filtered_features[position[1]].qualifiers["ApEinfo_revcolor"][0]

            if len(position) == 0:
                if level[0] == -1:
                    sequence_buffer = sequence_buffer + goi[i]
                else:
                    sequence_buffer = "%s<span style='font-weight:normal'><font color='black'>%s" % (
                        sequence_buffer, goi[i])
                    layer = layer + 1
                level[0:1] = [-1, -1]
            elif len(position) == 1:
                if level[1] != -1:  # Bugfix to: COLOR DOESN'T RETURN TO ORIGINAL WHEN 2nd LEVEL CLOSED
                    sequence_buffer = "%s</font><span style='font-weight:normal'>" % (
                        sequence_buffer)
                if position[0] == level[0]:
                    sequence_buffer = sequence_buffer + goi[i]
                else:
                    layer = layer + 1
                    level[0] = position[0]
                    sequence_buffer = "%s<span style='font-weight:normal'><font color='%s'>%s" % (
                        sequence_buffer, col[0], goi[i])
                    annotation_buffer = "%s<span style='font-weight:normal'><font color='%s'> %s" % (
                        annotation_buffer, col[0], filtered_features[position[0]].qualifiers["label"][0])
                level[1] = -1
                continue
            elif len(position) >= 2:
                if position[0] == level[0]:
                    if position[1] == level[1]:
                        sequence_buffer = sequence_buffer + goi[i]
                    else:
                        level[1] = position[1]
                        layer = layer + 1
                        sequence_buffer = "%s<span style='font-weight:bold'><font color='%s'>%s" % (
                            sequence_buffer, col[1], goi[i])
                        annotation_buffer = "%s<span style='font-weight:bold'><font color='%s'> %s" % (
                            annotation_buffer, col[1], filtered_features[position[1]].qualifiers["label"][0])
                else:
                    level[0], level[1] = position[0], position[1]
                    layer = layer + 1
                    sequence_buffer = "%s<span style='font-weight:normal'><font color='%s'>%s" % (
                        sequence_buffer, col[0], goi[i])
                    annotation_buffer = "%s<span style='font-weight:normal'><font color='%s'> %s" % (
                        annotation_buffer, col[0], filtered_features[position[0]].qualifiers["label"][0])
                    annotation_buffer = "%s<span style='font-weight:bold'><font color='%s'> %s" % (
                        annotation_buffer, col[1], filtered_features[position[1]].qualifiers["label"][0])
            else:  # len(position)>2
                sequence_buffer = sequence_buffer + goi[i]
        sequence_buffer = sequence_buffer + layer * "</font>"
        annotation_buffer = annotation_buffer + layer * "</font>"

        # Bugfix extends features to multiples of three, to color SNPs. In principle, this makes the checks below redundand
        for index in range(0, len(featurecomprehension), 3):
            buffer = []
            for comprehension in featurecomprehension[index:index + 2]:
                for element in comprehension:
                    if element not in buffer:
                        buffer.append(element)
            featurecomprehension[index] = buffer
            featurecomprehension[index + 1] = buffer
            featurecomprehension[index + 2] = buffer

        # labeling of protein sequence
        level = [-1, -1]
        layer = 0
        translation_buffer = ""
        protein_buffer = ""
        for i, position in enumerate(featurecomprehension):
            if not labelcolor:
                if len(position) > 0:
                    col[0] = colors[position[0]]
                    if len(position) > 1:
                        col[1] = colors[position[1]]
            else:
                if len(position) > 0:
                    col[0] = filtered_features[position[0]].qualifiers["ApEinfo_revcolor"][0]
                    if len(position) > 1:
                        col[1] = filtered_features[position[1]].qualifiers["ApEinfo_revcolor"][0]

            if len(translation_buffer) % 3 != 0:  # No changes in color before triplet filled. Should give a warning?
                translation_buffer = translation_buffer + goi[i]
            elif len(position) == 0:
                if level[0] == -1:
                    translation_buffer = translation_buffer + goi[i]
                else:
                    protein_buffer = "%s%s<span style='font-weight:normal'><font color='black'>" % (
                        protein_buffer, Seq(translation_buffer).translate())
                    translation_buffer = ""
                    translation_buffer = translation_buffer + goi[i]
                    layer = layer + 1
                level[0], level[1] = [-1, -1]
            elif len(position) == 1:
                if level[1] != -1:  # Bugfix to: COLOR DOESN'T RETURN TO ORIGINAL WHEN 2nd LEVEL CLOSED
                    protein_buffer = "%s%s<span style='font-weight:normal'><font color='%s'>" % (
                        protein_buffer, Seq(translation_buffer).translate(), col[0])
                    translation_buffer = ""
                    translation_buffer = translation_buffer + goi[i]
                    layer = layer + 1
                    level[0] = position[0]
                elif position[0] == level[0]:
                    translation_buffer = translation_buffer + goi[i]
                else:
                    protein_buffer = "%s%s<span style='font-weight:normal'><font color='%s'>" % (
                        protein_buffer, Seq(translation_buffer).translate(), col[0])
                    translation_buffer = ""
                    translation_buffer = translation_buffer + goi[i]
                    layer = layer + 1
                    level[0] = position[0]
                level[1] = -1
            elif len(position) >= 2:
                if position[0] == level[0]:
                    if position[1] == level[1]:
                        translation_buffer = translation_buffer + goi[i]
                    else:
                        protein_buffer = "%s%s<span style='font-weight:bold'><font color='%s'>" % (
                            protein_buffer, Seq(translation_buffer).translate(), col[1])
                        translation_buffer = ""
                        translation_buffer = translation_buffer + goi[i]
                        layer = layer + 1
                        level[1] = position[1]
                else:
                    level[0], level[1] = position[0], position[1]
                    layer = layer + 1
                    protein_buffer = "%s%s<span style='font-weight:normal'><font color='%s'>" % (
                        protein_buffer, Seq(translation_buffer).translate(), col[0])
                    protein_buffer = "%s%s<span style='font-weight:bold'><font color='%s'>" % (
                        protein_buffer, Seq(translation_buffer).translate(), col[1])
                    translation_buffer = ""
                    translation_buffer = translation_buffer + goi[i]
            else:  # len(position)>2
                translation_buffer = translation_buffer + goi[i]
        protein_buffer = protein_buffer + str(Seq(translation_buffer).translate()) + layer * "</font>"
        protein_buffer = protein_buffer.replace("*", "-")

    else:
        sequence_buffer = str(goi)
        protein_buffer = str(goi.translate())
        annotation_buffer = ""

    # Resolving template to result file
    templatefile = os.path.join(os.getcwd(), "seq_align", "templates", "seq_val_raw_template.htm")
    resultfile = os.path.join(os.getcwd(), html_dict["construct"] + html_dict["plasmid"] + ".htm")
    with open(templatefile, "r") as fh:
        with open(resultfile, "w+") as rh:
            for line in fh:
                split, join = line.split("#"), ""
                for part in split:  # you could just do this with {key}.format(dict)
                    if part == "seqtemplate":
                        # In code as this has to resolve multiple times. Could be own function
                        seqbuffer = ""
                        for element in alignment_list:
                            ab1 = element["ab1"]
                            alignment = element["alignment"]
                            traces = element["traces"]
                            seqformat = seqtemplate
                            seqformat = seqformat.replace("#SEQ_NAME#", ab1.name)
                            seqformat = seqformat.replace("#SEQ_PRIMER#", ab1.annotations["abif_raw"]["CMNT1"])
                            seqformat = seqformat.replace("#SEQ_DATE#", ab1.annotations["abif_raw"]["RUND1"])
                            seqformat = seqformat.replace("#SEQ_READ#", str(ab1.seq))
                            try:
                                alignment_range = "%s-%s, %s-%s" % (
                                    alignment["que_range"][0], traces[0]["position"], traces[-1]["position"],
                                    alignment["que_range"][1])
                                if alignment["que_range"][1] == len(goi):
                                    alignment_range = alignment_range + " (end)"
                            except:
                                alignment_range = "%s-%s" % (alignment["que_range"][0], alignment["que_range"][1])
                                if alignment["que_range"][1] == len(goi):
                                    alignment_range = alignment_range + " (end)"
                            seqformat = seqformat.replace("#BLAST_VAL#", alignment_range)
                            seqformat = seqformat.replace("#BLAST#", alignment["html_tiny_blast"])
                            for trace in traces:
                                seqformat = seqformat.replace("#TRACE#",
                                                              "<img src='%s' width='300'> #TRACE#" % trace[
                                                                  "figurefile"])
                                seqformat = seqformat.replace("#TRACE_REG#", trace["mutation"] + " #TRACE_REG#")
                            if len(traces) == 0:
                                seqformat = seqformat.replace("#TRACE_REG#", "None")
                            else:
                                seqformat = seqformat.replace("#TRACE_REG#", "")
                            seqformat = seqformat.replace("#TRACE#", "")
                            seqbuffer = seqbuffer + ab1.annotations["abif_raw"][
                                "CMNT1"] + " gives " + alignment_range + "<br>&emsp;"
                            join = join + seqformat
                    elif part == "SEQSUM":
                        join = join + seqbuffer
                    elif part == "sequence":
                        join = join + sequence_buffer
                    elif part == "PROTEIN":
                        join = join + protein_buffer
                    elif part == "ANNOTATION":
                        join = join + annotation_buffer
                    elif part in html_dict:
                        if isinstance(html_dict[part], list):
                            join = join + str(html_dict[part][0])  # fix this later
                        else:
                            join = join + str(html_dict[part])
                    else:
                        join = join + part
                rh.write(join)
    print("Task processed. Output %s \n" % resultfile)

if __name__ == '__main__':
    username = "Your name"  # Put your name here
    verbose = True  # Use this if you want extra output information
    labelcolor = False  # Use this if you want the labels to be colored according to genbank file
    colors = ["#e6194b", "#3cb44b", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
              "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080",
              "#808080"]  # colors = ["blue", "green", "red", "purple", "darkorange", "mustard", "teal", "brown", "cyan", "magenta",]  # redefine color scheme
    # wd = os.getcwd()  # Use this to set directory. By default, terminal directory. Change e.g.
    wd = "C:\\Users\\Arne Scheu\\Desktop\\qol"
    master()
