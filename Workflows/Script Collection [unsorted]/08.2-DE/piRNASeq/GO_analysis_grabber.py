import os
import pandas as pd
import json
import datetime
import sys

def define_up_down_range(log2foldchange: float):
    change = abs(log2foldchange * 10)
    entry = "<.1"
    for x in range(1, 9):
        if x < change < x + 1:
            entry = f".{x}-.{x + 1}"
    if 9 < change:
        entry = ".9+"
    return entry

# iterate over the names instead of doing two for loops...
filenames = ["polyIC", "ecoli"]

cond_dict, ecoli_all, polyIC_all = dict(), dict(), dict()

with open("/Work2/Planarian_genome_2021/SMEST2G_2GO/SMEST_to_SMESG.tsv") as tsv:
    smesg2smest = dict()
    for line in tsv:
        smesg = line.split()[1]
        smest = line.split()[0]
        if smesg not in smesg2smest:
            smesg2smest[smesg] = []
        smesg2smest[smesg].append(smest)

if __name__ == "__main__":
    try:
        folder = sys.argv[1]
    except:
        print("No folder specified... Folder needed for SMESG annotation")
        exit()

    cluster2smesg=json.load(open("pdens0.3_swsize3000_clsize3000.cluster2smesg.json", "r"))
    cluster2smesg=cluster2smesg.get(folder)
    os.chdir(folder)

    try:
        print("Creating folder Analysis")
        os.mkdir("./Analysis")
    except:
        print("Folder Analysis already exists.")

    # Taking transcript, log2FoldChange, pvalue from the table and converting into a dict()

    for typ_cond in filenames:
        df = pd.read_table(open(f"de_unique_{typ_cond}.txt", "r"), names=["transcript", "baseMean", "log2FoldChange",
                                                                          "lfcSE", "stat", "pvalue", "padj"],
                           skiprows=1)
        df = df[["transcript", "log2FoldChange", "pvalue"]]
        transcripts_only = {
            "UP": {"<.1": {}, ".1-.2": {}, ".2-.3": {}, ".3-.4": {}, ".4-.5": {}, ".5-.6": {}, ".6-.7": {}, ".7-.8": {},
                   ".8-.9": {},
                   ".9+": {}},
            "DOWN": {"<.1": {}, ".1-.2": {}, ".2-.3": {}, ".3-.4": {}, ".4-.5": {}, ".5-.6": {}, ".6-.7": {},
                     ".7-.8": {}, ".8-.9": {},
                     ".9+": {}}}
        upreg = dict()
        downreg = dict()
        for index, transcript in df.iterrows():
            transcript_id = transcript["transcript"]
            transcript_log2foldchange = transcript["log2FoldChange"]
            transcript_pValue = transcript["pvalue"]
            smests = [smesg2smest.get(x) for x in cluster2smesg.get(transcript_id)]
            locals()[f"{typ_cond}_all"][transcript_id] = {"log2FoldChange": transcript_log2foldchange,
                                                          "pvalue": transcript_pValue,
                                                          "SMESGs": cluster2smesg.get(transcript_id),
                                                          "SMESTs":smests}
            updown = "DOWN"
            if transcript_log2foldchange > 0:
                updown = "UP"
            entry_range = define_up_down_range(transcript_log2foldchange)
            transcripts_only[updown][entry_range][transcript_id] = {"log2FoldChange": transcript_log2foldchange,
                                                                    "pvalue": transcript_pValue,
                                                                    "SMESGs": cluster2smesg.get(transcript_id),
                                                                    "SMESTs":smests}

        cond_dict[typ_cond] = transcripts_only

    # because I love jsons... storing the information into a json file.
    # {Ecoli: {Transcript : {log2FoldChange, pvalue}, ...}, PolyIC {Transcript : {log2FoldChange, pvalue}, ...}}
    json.dump(cond_dict, open("./Analysis/de_unique_transcriptsonly.json", "w"), indent=4)

    polyIC_vs_ecoli_sharedtranscripts = dict()
    ecoli_specific, polyIC_specific = [{
        "UP": {"<.1": {}, ".1-.2": {}, ".2-.3": {}, ".3-.4": {}, ".4-.5": {}, ".5-.6": {}, ".6-.7": {}, ".7-.8": {},
               ".8-.9": {},
               ".9+": {}},
        "DOWN": {"<.1": {}, ".1-.2": {}, ".2-.3": {}, ".3-.4": {}, ".4-.5": {}, ".5-.6": {}, ".6-.7": {}, ".7-.8": {},
                 ".8-.9": {},
                 ".9+": {}}
    }
        for x in range(0, 2)]

    for transcript, info in ecoli_all.items():
        ecoli_log2 = info.get("log2FoldChange")
        ecoli_pvalue = info.get("pvalue")
        polyIC_updown, ecoli_updown = "DOWN", "DOWN"
        smests = [smesg2smest.get(x) for x in cluster2smesg.get(transcript)]
        if ecoli_log2 > 0:
            ecoli_updown = "UP"
        if transcript in polyIC_all:
            polyIC_log2 = polyIC_all.get(transcript).get("log2FoldChange")
            polyIC_pvalue = polyIC_all.get(transcript).get("pvalue")
            if polyIC_log2 > 0:
                polyIC_updown = "UP"
            vs_state = f"{ecoli_updown}:{polyIC_updown}"
            polyIC_vs_ecoli_sharedtranscripts[transcript] = {"ecoli:polyIC": vs_state,
                                                             "ecoli": {"state": ecoli_updown,
                                                                       "log2FoldChange": ecoli_log2,
                                                                       "pvalue": ecoli_pvalue,
                                                          "SMESGs": cluster2smesg.get(transcript),
                                                                        "SMESTs":smests},
                                                             "polyIC": {"state": polyIC_updown,
                                                                        "log2FoldChange": polyIC_log2,
                                                                        "pvalue": polyIC_pvalue,
                                                          "SMESGs": cluster2smesg.get(transcript),
                                                                        "SMESTs":smests}}

        else:
            inbetween = define_up_down_range(ecoli_log2)
            ecoli_specific[ecoli_updown][inbetween][transcript] = {"log2FoldChange": ecoli_log2, "pvalue": ecoli_pvalue,
                                                          "SMESGs": cluster2smesg.get(transcript),"SMESTs":smests}

    for transcript, info in polyIC_all.items():
        smests = [smesg2smest.get(x) for x in cluster2smesg.get(transcript)]
        polyIC_log2 = info.get("log2FoldChange")
        polyIC_pvalue = info.get("pvalue")
        polyIC_updown = "DOWN"
        if polyIC_log2 > 0:
            polyIC_updown = "UP"
        if transcript not in ecoli_all:
            inbetween = define_up_down_range(polyIC_log2)
            polyIC_specific[polyIC_updown][inbetween][transcript] = {"log2FoldChange": polyIC_log2,
                                                                     "pvalue": polyIC_pvalue,
                                                          "SMESGs": cluster2smesg.get(transcript),"SMESTs":smests}

    ## Check, if any are upregulated vs downregulated depending on the condition
    for transcript, information in polyIC_vs_ecoli_sharedtranscripts.items():
        vs_state = information.get("ecoli:polyIC")
        # if vs_state not in ["UP:UP", "DOWN:DOWN"]:
        #    print(transcript) # We conclude,... if it's up/downregulated, its the same in both conditions.......

    json.dump(polyIC_vs_ecoli_sharedtranscripts, open("Analysis/de_shared_transcriptsonly.json", "w"), indent=4)
    json.dump(polyIC_specific, open("Analysis/de_polyIC_specific_transcriptsonly.json", "w"), indent=4)
    json.dump(ecoli_specific, open("Analysis/de_ecoli_specific_transcriptsonly.json", "w"), indent=4)

    with open("Analysis/de_unique_statistics.txt", "w") as handle:
        handle.write(f"{datetime.datetime.now().strftime('%a %H:%M - %Y-%m-%d')}\n"
                     f"------------------------------\n"
                     f"Ecoli Total: {len(ecoli_all)}\n"
                     f"Unique to Ecoli (not in PolyIC): {len(ecoli_all) - len(polyIC_vs_ecoli_sharedtranscripts)}\n\n"
                     f"PolyIC Total: {len(polyIC_all)}\n"
                     f"Unique to PolyIC (not in Ecoli): {len(polyIC_all) - len(polyIC_vs_ecoli_sharedtranscripts)}\n\n"
                     f"Shared: {len(polyIC_vs_ecoli_sharedtranscripts)}")


    all_unique = json.load(open("Analysis/de_unique_transcriptsonly.json", "r"))

    polyIC_SMESTs_UP = list()
    polyIC_SMESTs_DOWN = list()
    for up_down, range_dict in all_unique.get("polyIC").items():
        for range in range_dict.items():
            for cluster, information in range[1].items():
                for smests in information.get("SMESTs"):
                    if smests not in locals()[f"polyIC_SMESTs_{up_down}"]:
                        locals()[f"polyIC_SMESTs_{up_down}"] += smests
    with open(f"Analysis/{folder}_unique_transcript_IDs_polyIC_UP.txt", "w") as out_up, \
        open(f"Analysis/{folder}_unique_transcript_IDs_polyIC_DOWN.txt", "w") as out_down, \
        open(f"Analysis/{folder}_unique_transcript_IDs_polyIC_all.txt", "w") as out_all:
        for n in sorted(polyIC_SMESTs_UP):
            out_up.write(n + "\n")
            out_all.write(n + "\n")
        for n in sorted(polyIC_SMESTs_DOWN):
            out_down.write(n + "\n")
            out_all.write(n + "\n")

    ecoli_SMESTs_UP = list()
    ecoli_SMESTs_DOWN = list()
    for up_down, range_dict in all_unique.get("ecoli").items():
        for range in range_dict.items():
            for cluster, information in range[1].items():
                for smests in information.get("SMESTs"):
                    if smests not in locals()[f"ecoli_SMESTs_{up_down}"]:
                        locals()[f"ecoli_SMESTs_{up_down}"] += smests
    with open(f"Analysis/{folder}_unique_transcript_IDs_ecoli_UP.txt", "w") as out_up, \
        open(f"Analysis/{folder}_unique_transcript_IDs_ecoli_DOWN.txt", "w") as out_down, \
        open(f"Analysis/{folder}_unique_transcript_IDs_ecoli_all.txt", "w") as out_all:
        for n in sorted(ecoli_SMESTs_UP):
            out_up.write(n + "\n")
            out_all.write(n + "\n")
        for n in sorted(ecoli_SMESTs_DOWN):
            out_down.write(n + "\n")
            out_all.write(n + "\n")