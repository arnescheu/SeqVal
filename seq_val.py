import csv
import datetime
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib.lines import Line2D


def master():
    # Defining which files to handle from task.txt file
    task_list = define_task(os.path.join(wd, "task.txt"))

    for task in task_list:
        # Defining gene of interest based on genbank and region or sequence field.
        # Returns a dictionary with sequence, start, end (+1 to array index)
        region = define_goi(task)
        goi, start, end = region["goi"], region["start"], region["end"]  # TODO: Make less ugly return

        # define the sequence traces belonging to gb file
        ab1_list = []
        for key in task:
            if "sequencing" in key:
                if task[key] is not None:
                    split = task[key].split(";")
                    for path in split:
                        ab1_list.append(path)
        # defaulting to all ab1 in directory if none defined
        if not ab1_list:
            print("No .ab1 traces defined for %s. Default to all traces." % task["genbank"])
            for files in os.walk(wd):
                if ".ab1" in files:
                    ab1_list.append(files)

        # For each sequence trace, creating corresponding alignment in nblast
        # Then creates image of traces at mismatch positions
        alignment_list = []
        for ab1_path in ab1_list:
            if local_nblast(goi, ab1_path) != 0:
                continue
            alignment = interpret_blast()
            ab1 = SeqIO.read(ab1_path, 'abi')
            traces = handle_alignment(ab1, alignment)
            alignment_list.append({"ab1": ab1, "alignment": alignment, "traces": traces})

        # Create expasy, removing first residue ([3:]
        with  open(os.path.join(wd, "seq_align", "temp", "nomet_expasy.html"), "w+") as ex:
            ex.write(expasy(goi[3:].translate(to_stop=True)).decode("utf-8"))

        # Compiles all previous data
        task["user"] = username
        if "record" in region:
            html_love(task, region["goi"], alignment_list, region["record"])
        else:
            html_love(task, region["goi"], alignment_list)


def define_task(taskfile):
    task_list = []
    with open(taskfile, "r") as csvfile:  # put in header, then to disct with all the necessary elements
        task_reader = csv.DictReader(csvfile, delimiter="\t", quotechar="#")
        # https://stackoverflow.com/questions/15865750/dictreader-change-keys-to-upper
        for row in task_reader:
            task_list.append(row)
    return task_list


def define_goi(task):  # not bothering with -1 strand
    if not task["genbank"] is None:
        try:
            record = SeqIO.read(task["genbank"], "genbank")
            for feature in record.features:
                if "GOI" in feature.qualifiers["label"]:  # define label name by genebank file?
                    start, end = int(feature.location.start), int(feature.location.end)
                    return {"record": record, "goi": record.seq[start - 1:end - 1], "start": start, "end": end}
            if not task["region"] is None:  # define sequence from region
                start = int(task["region"].split(":")[0])
                end = int(task["region"].split(":")[1])
                return {"record": record, "goi": record.seq[start - 1:end - 1], "start": start, "end": end}
            elif not task["sequence"] is None:  # define region from sequence
                goi = task["sequence"]
                start = int(record.seq.find(goi)) - 1
                return {"record": record, "goi": goi, "start": start, "end": start + len(goi)}
            else:
                print("GOI not defined by region/sequence/annotation. Default to genbank region")
                return {"record": record, "goi": str(record.seq), "start": 1, "end": len(record.seq)}
        except FileNotFoundError:
            print(task["genebank"] + " not found.")
            if not task["sequence"] is None:
                return {"goi": task["sequence"], "start": 1, "end": len(task["sequence"])}
            else:
                return False
    elif not task["sequence"] is None:
        return {"goi": task["sequence"], "start": 1, "end": len(task["sequence"])}
    else:
        print("Task lacking either genbank with region and sequence.")
        return False


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
            'blastn -subject "%s" -query %s -dust no -outfmt 3 >%s' % (
        seq_temp, query_temp, out_temp))  # output 7 creates xml biopython blast output. Biopython can read
    return os.system(command)


def interpret_blast(
        blastfile="blast_temp.txt"):  # converts blast output so subject/query. Very bad practice, could just use other modules
    with open(os.path.join(os.getcwd(), "seq_align", "temp", blastfile), "r") as raw_blast:
        start, space = False, 0
        subject = {"query": "", "sequence": "", "misalign": []}
        subject["tiny_blast"] = ""
        subject["seq_range"] = [False, False]
        for line in raw_blast:
            space += 1
            if "Query_1" in line:
                start, space = True, 0
                subject["query"] = subject["query"] + list(filter(None, line.strip("\n").split(" ")))[2]
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
                    subject["tiny_blast"] = subject["tiny_blast"] + line.replace(linesplit[2], linesplit_rec)

                elif space > 2:  # stop after last query:subject block
                    break
                else:
                    subject["tiny_blast"] = subject["tiny_blast"] + line

        for i, character in enumerate(subject["sequence"]):
            if character != ".":
                subject["misalign"].append(i)
        style = "<font color='red'>"
        subject["html_tiny_blast"] = subject["tiny_blast"].replace("\n", "<br>").replace(" ", "&nbsp;").replace(
                "#style#",
                style)
    return subject


def handle_alignment(ab1, alignment):
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

    """ This produces the same guide, using a different method
    guide=list("."*len(alignment["sequence"]))
    i=0
    for el in alignment["sequence"].split("-"):
        i = i + len(el)
        if i != len(guide):
            guide[i] = "-"
        i = i + 1
    i=0
    for el in alignment["query"].split("-"):
        i=i+len(el)
        if i != len(guide):
            guide[i]="+"
        i=i+1
    guide="".join(guide)
    #print("guide_1",len(guide),guide)
    """

    offset = 0
    for relative_misalign in alignment["misalign"]:
        if alignment["seq_range"][0] < alignment["seq_range"][1]:
            alignment["orientation"] = 1
            traces.append(chromatogram(relative_misalign, offset, ab1, alignment))
            if guide[relative_misalign] == ".":
                pass
            elif guide[relative_misalign] == "+":
                offset = offset + 1
            elif guide[relative_misalign] == "-":
                offset = offset - 1
        else:
            alignment["orientation"] = -1
            if guide[relative_misalign] == ".":
                pass
            elif guide[relative_misalign] == "+":
                offset = offset + 1
            elif guide[relative_misalign] == "-":
                offset = offset - 1
            traces.append(chromatogram(relative_misalign, offset, ab1, alignment))
    return traces


def chromatogram(relative_misalign, offset, ab1, alignment):
    if alignment["orientation"] == 1:
        rel = relative_misalign + offset - 1
        absolute_position = alignment["seq_range"][0] + rel
    else:
        rel = -relative_misalign - offset
        absolute_position = alignment["seq_range"][0] + rel

    readposition = ab1.annotations["abif_raw"]["PLOC1"][absolute_position]

    print("actual position and read", str(absolute_position), ab1.seq[absolute_position - 2:absolute_position + 2],
          ab1.seq[absolute_position])  # 1227,1224
    print("relative_misalign and query, subject", relative_misalign, alignment["query"][relative_misalign],
          alignment["sequence"][relative_misalign])
    print("read position", readposition)

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
        print("Index end of chromatogram, defaulting to maximum")
        xmax = len(ab1.annotations["abif_raw"]["DATA9"])
    plt.xlim(xmin, xmax)
    plt.axvline(x=readposition)
    plt.ylim(0, ymax)

    ax1 = plt.gca()
    ax1.set_xlabel('Chromatogram')
    ax = plt.subplot()

    plt.title("Misaligned " + alignment["query"][relative_misalign] + str(relative_misalign) + alignment["sequence"][
        relative_misalign] + ", ab1 " + str(
            readposition) + " offset " + str(offset))  # ,  y=1.125)
    # plt.title(alignment["query"][relative_misalign-10:relative_misalign+10],y=1.125)
    # TODO alignment Query:Sequence

    # ax2=ax1.twiny()
    # ax2.set_xlabel('Query nt')

    # TODO fix pointer location
    # Adding pointers on 2nd axis
    visible_ticks = []
    for r in range(absolute_position - 8, absolute_position + 10, 1):
        visible_ticks.append(ab1.annotations["abif_raw"]["PLOC1"][r] - xmin)
    print(visible_ticks)
    # ax2.set_xticks(visible_ticks)
    # ax2.set_xticklabels(range(relative_misalign-8,relative_misalign+10,1))
    # ax2.plot(ab1.annotations["abif_raw"]["PLOC1"])

    # Adding legend for G/A/T/C
    # stolen from somewhere on stackoverflow
    custom_lines = [Line2D([0], [0], color=clr["DATA9"], lw=4),
                    Line2D([0], [0], color=clr["DATA10"], lw=4),
                    Line2D([0], [0], color=clr["DATA11"], lw=4),
                    Line2D([0], [0], color=clr["DATA12"], lw=4), ]
    ax.legend(custom_lines, ["G", "A", "T", "C"])

    # saving file
    figname = os.path.join("seq_align", "results", "traces",
                           alignment["query"][relative_misalign] + str(relative_misalign) + alignment["sequence"][
                               relative_misalign] + "_" + ab1.name + ".png")
    plt.savefig(os.path.join(os.getcwd(), figname))
    plt.show()

    return [figname,
            alignment["query"][relative_misalign] + str(relative_misalign) + alignment["sequence"][relative_misalign]]


def expasy(sequence):  # remove initial methionine
    url = "https://web.expasy.org/cgi-bin/protparam/protparam"
    payload = {'sequence': sequence}
    r = requests.get(url, params=payload)
    return r.content


def html_love(html_dict, goi, alignment_list, record=None):
    html_dict["date_seqval"] = "%s-%s-%s" % (datetime.datetime.now().year, str(datetime.datetime.now().month).zfill(2),
                                             str(datetime.datetime.now().day).zfill(2))
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
    sorted_features = sorted(record.features, key=lambda x: x.location.start)
    with open(os.path.join(os.getcwd(), "seq_align", "templates", "seqtemplate.txt"), "r") as fh:
        seqtemplate = fh.read()
    colors = ["blue", "green", "purple", "red", "orange", "pink", "teal", "brown", "cyan", "magenta"]

    featurecomprehension = []
    for i in range(0, len(goi)):
        featurecomprehension.append([])
    try:
        start = int(html_dict["region"].split(":")[0])
    except KeyError:
        start = 0
    filtered_features = []
    i = 0
    for feature in sorted_features:
        if not (feature.type == "CDS" or feature.type == "Protein"):
            continue
        if feature.location.strand == -1:
            continue
        featurestart = int(feature.location.start) - start + 1
        featureend = int(feature.location.end) - start + 1
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
    for i, position in enumerate(featurecomprehension):
        if len(position) == 0:
            if level[0] == -1:
                sequence_buffer = sequence_buffer + goi[i]
            else:
                sequence_buffer = "%s<span style='font-weight:normal'><font color='black'>%s" % (
                    sequence_buffer, goi[i])
                layer = layer + 1
            level[0:1] = [-1, -1]
        elif len(position) == 1:
            if position[0] == level[0]:
                sequence_buffer = sequence_buffer + goi[i]
            else:
                layer = layer + 1
                level[0] = position[0]
                sequence_buffer = "%s<span style='font-weight:normal'><font color='%s'>%s" % (
                    sequence_buffer, colors[position[0]], goi[i])
                annotation_buffer = "%s<span style='font-weight:normal'><font color='%s'> %s" % (
                    annotation_buffer, colors[position[0]], filtered_features[position[0]].qualifiers["label"][0])
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
                        sequence_buffer, colors[position[1]], goi[i])
                    annotation_buffer = "%s<span style='font-weight:bold'><font color='%s'> %s" % (
                        annotation_buffer, colors[position[1]], filtered_features[position[1]].qualifiers["label"][0])
            else:
                level[0], level[1] = position[0], position[1]
                layer = layer + 1
                sequence_buffer = "%s<span style='font-weight:bold'><font color='%s'>%s" % (
                    sequence_buffer, colors[position[0]], goi[i])
                annotation_buffer = "%s<span style='font-weight:normal'><font color='%s'> %s" % (
                    annotation_buffer, colors[position[0]], filtered_features[position[0]].qualifiers["label"][0])
                annotation_buffer = "%s<span style='font-weight:bold'><font color='%s'> %s" % (
                    annotation_buffer, colors[position[1]], filtered_features[position[1]].qualifiers["label"][0])
        else:  # len(position)>2
            sequence_buffer = sequence_buffer + goi[i]
    sequence_buffer = sequence_buffer + layer * "</font>"
    annotation_buffer = annotation_buffer + layer * "</font>"

    # labeling of protein sequence
    level = [-1, -1]
    layer = 0
    translation_buffer = ""
    protein_buffer = ""
    for i, position in enumerate(featurecomprehension):
        if len(translation_buffer) % 3 != 0:
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
            if position[0] == level[0]:
                translation_buffer = translation_buffer + goi[i]
            else:
                protein_buffer = "%s%s<span style='font-weight:normal'><font color='%s'>" % (
                    protein_buffer, Seq(translation_buffer).translate(), colors[position[0]])
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
                        protein_buffer, Seq(translation_buffer).translate(), colors[position[1]])
                    translation_buffer = ""
                    translation_buffer = translation_buffer + goi[i]
                    layer = layer + 1
                    level[1] = position[1]
            else:
                level[0], level[1] = position[0], position[1]
                layer = layer + 1
                protein_buffer = "%s%s<span style='font-weight:normal'><font color='%s'>" % (
                    protein_buffer, Seq(translation_buffer).translate(), colors[position[0]])
                protein_buffer = "%s%s<span style='font-weight:bold'><font color='%s'>" % (
                    protein_buffer, Seq(translation_buffer).translate(), colors[position[1]])
                translation_buffer = ""
                translation_buffer = translation_buffer + goi[i]
        else:  # len(position)>2
            translation_buffer = translation_buffer + goi[i]
    protein_buffer = protein_buffer + str(Seq(translation_buffer).translate()) + layer * "</font>"
    protein_buffer = protein_buffer.replace("*", "-")

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
                            alignment_range = "1-%s,%s-end" % (min(alignment["misalign"]), max(alignment["misalign"]))
                            seqformat = seqformat.replace("#BLAST_VAL#", alignment_range)
                            seqformat = seqformat.replace("#BLAST#", alignment["html_tiny_blast"])
                            for trace in traces:
                                seqformat = seqformat.replace("#TRACE#",
                                                              "<img src='%s' width='300'> #TRACE#" % trace[0])
                                seqformat = seqformat.replace("#TRACE_REG#", trace[1] + " #TRACE_REG#")
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


if __name__ == '__main__':
    wd = "C:\\Users\\Arne Scheu\\Desktop\\qol"
    username = "Arne Scheu"
    os.chdir(wd)
    master()
