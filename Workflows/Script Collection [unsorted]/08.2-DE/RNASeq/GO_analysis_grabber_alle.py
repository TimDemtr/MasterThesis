import os
import pandas as pd
import json
import datetime
import sys

def define_up_down_range(log2foldchange: float):
    change = abs(log2foldchange)
    entry = "<1"
    for x in range(1, 9):
        if x < change < x + 1:
            entry = f"{x}-{x + 1}"
    if 9 < change:
        entry = "9+"
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
        cutoff = (sys.argv[1])
    except:
        cutoff = 0

    print(f"log2foldchange cutoff: {cutoff}")

    try:
        print("Creating folder Analysis_alle")
        os.mkdir("./Analysis_alle")
    except:
        print("Folder Analysis already exists.")

    # Taking transcript, log2FoldChange, pvalue from the table and converting into a dict()

    for typ_cond in filenames:
        df = pd.read_table(open(f"de_alle_{typ_cond}.txt", "r"), names=["transcript", "baseMean", "log2FoldChange",
                                                                          "lfcSE", "stat", "pvalue", "padj"],
                           skiprows=1)
        df = df[["transcript", "log2FoldChange", "pvalue"]]
        transcripts_only = {
            "UP": {"<1": {}, "1-2": {}, "2-3": {}, "3-4": {}, "4-5": {}, "5-6": {}, "6-7": {}, "7-8": {},
                   "8-9": {},
                   "9+": {}},
            "DOWN": {"<1": {}, "1-2": {}, "2-3": {}, "3-4": {}, "4-5": {}, "5-6": {}, "6-7": {},
                     "7-8": {}, "8-9": {},
                     "9+": {}}}
        upreg = dict()
        downreg = dict()
        for index, transcript in df.iterrows():
            transcript_id = transcript["transcript"]
            transcript_log2foldchange = transcript["log2FoldChange"]
            if transcript_log2foldchange == "NaN":
                continue
            transcript_pValue = transcript["pvalue"]
            locals()[f"{typ_cond}_all"][transcript_id] = {"log2FoldChange": transcript_log2foldchange,
                                                          "pvalue": transcript_pValue,
                                                                        "SMESTs":smesg2smest.get(transcript_id)}
            updown = "DOWN"
            if transcript_log2foldchange > 0:
                updown = "UP"
            entry_range = define_up_down_range(transcript_log2foldchange)
            transcripts_only[updown][entry_range][transcript_id] = {"log2FoldChange": transcript_log2foldchange,
                                                                    "pvalue": transcript_pValue,
                                                                        "SMESTs":smesg2smest.get(transcript_id)}

        cond_dict[typ_cond] = transcripts_only

    # because I love jsons... storing the information into a json file.
    # {Ecoli: {Transcript : {log2FoldChange, pvalue}, ...}, PolyIC {Transcript : {log2FoldChange, pvalue}, ...}}
    json.dump(cond_dict, open("./Analysis_alle/de_alle_transcriptsonly.json", "w"), indent=4)

    polyIC_vs_ecoli_sharedtranscripts = dict()
    ecoli_specific, polyIC_specific = [{
        "UP": {"<1": {}, "1-2": {}, "2-3": {}, "3-4": {}, "4-5": {}, "5-6": {}, "6-7": {},
                     "7-8": {}, "8-9": {},
                     "9+": {}},
        "DOWN": {"<1": {}, "1-2": {}, "2-3": {}, "3-4": {}, "4-5": {}, "5-6": {}, "6-7": {},
                     "7-8": {}, "8-9": {},
                     "9+": {}}
    }
        for x in range(0, 2)]

    for transcript, info in ecoli_all.items():
        ecoli_log2 = info.get("log2FoldChange")
        ecoli_pvalue = info.get("pvalue")
        polyIC_updown, ecoli_updown = "DOWN", "DOWN"
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
                                                                       "SMESTs": smesg2smest.get(transcript)},
                                                             "polyIC": {"state": polyIC_updown,
                                                                        "log2FoldChange": polyIC_log2,
                                                                        "pvalue": polyIC_pvalue,
                                                                        "SMESTs":smesg2smest.get(transcript)}}

        else:
            inbetween = define_up_down_range(ecoli_log2)
            ecoli_specific[ecoli_updown][inbetween][transcript] = {"log2FoldChange": ecoli_log2, "pvalue": ecoli_pvalue,"SMESTs": smesg2smest.get(transcript)}

    for transcript, info in polyIC_all.items():
        polyIC_log2 = info.get("log2FoldChange")
        polyIC_pvalue = info.get("pvalue")
        polyIC_updown = "DOWN"
        if polyIC_log2 > 0:
            polyIC_updown = "UP"
        if transcript not in ecoli_all:
            inbetween = define_up_down_range(polyIC_log2)
            polyIC_specific[polyIC_updown][inbetween][transcript] = {"log2FoldChange": polyIC_log2,
                                                                     "pvalue": polyIC_pvalue,
                                                                     "SMESTs": smesg2smest.get(transcript)}

    ## Check, if any are upregulated vs downregulated depending on the condition
    for transcript, information in polyIC_vs_ecoli_sharedtranscripts.items():
        vs_state = information.get("ecoli:polyIC")
        # if vs_state not in ["UP:UP", "DOWN:DOWN"]:
        #    print(transcript) # We conclude,... if it's up/downregulated, its the same in both conditions.......

    json.dump(polyIC_vs_ecoli_sharedtranscripts, open("Analysis_alle/de_shared_transcriptsonly.json", "w"), indent=4)
    json.dump(polyIC_specific, open("Analysis_alle/de_polyIC_specific_transcriptsonly.json", "w"), indent=4)
    json.dump(ecoli_specific, open("Analysis_alle/de_ecoli_specific_transcriptsonly.json", "w"), indent=4)

    with open("Analysis_alle/de_alle_statistics.txt", "w") as handle:
        handle.write(f"{datetime.datetime.now().strftime('%a %H:%M - %Y-%m-%d')}\n"
                     f"------------------------------\n"
                     f"Ecoli Total: {len(ecoli_all)}\n"
                     f"Unique to Ecoli (not in PolyIC): {len(ecoli_all) - len(polyIC_vs_ecoli_sharedtranscripts)}\n\n"
                     f"PolyIC Total: {len(polyIC_all)}\n"
                     f"Unique to PolyIC (not in Ecoli): {len(polyIC_all) - len(polyIC_vs_ecoli_sharedtranscripts)}\n\n"
                     f"Shared: {len(polyIC_vs_ecoli_sharedtranscripts)}")


    all_unique = json.load(open("Analysis_alle/de_alle_transcriptsonly.json", "r"))
    try:
        os.mkdir(f"Analysis_alle/transcripts_only/")
        os.mkdir(f"Analysis_alle/transcripts_only/{cutoff}/")
    except:
        try:
            os.mkdir(f"Analysis_alle/transcripts_only/{cutoff}/")
        except:
            pass
    polyIC_SMESTs_UP = list()
    polyIC_SMESTs_DOWN = list()
    for up_down, range_dict in all_unique.get("polyIC").items():
        for range, smesgs in range_dict.items():
            for smesg,information in smesgs.items():
                if abs(information.get("log2FoldChange")) <= abs(float(cutoff)):
                    continue
                for smests in information.get("SMESTs"):
                    if smests not in locals()[f"polyIC_SMESTs_{up_down}"] and smests.endswith("1.1"):
                        locals()[f"polyIC_SMESTs_{up_down}"].append(smests)
    with open(f"Analysis_alle/transcripts_only/{cutoff}/alle_transcript_IDs_polyIC_UP_cutoff{cutoff}.txt", "w") as out_up, \
        open(f"Analysis_alle/transcripts_only/{cutoff}/alle_transcript_IDs_polyIC_DOWN_cutoff{cutoff}.txt", "w") as out_down, \
        open(f"Analysis_alle/transcripts_only/{cutoff}/alle_transcript_IDs_polyIC_all_cutoff{cutoff}.txt", "w") as out_all:
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
                if abs(information.get("log2FoldChange")) <= abs(float(cutoff)):
                    continue
                for smests in information.get("SMESTs"):
                    if smests not in locals()[f"ecoli_SMESTs_{up_down}"] and smests.endswith("1.1"):
                        locals()[f"ecoli_SMESTs_{up_down}"].append(smests)
    with open(f"Analysis_alle/transcripts_only/{cutoff}/alle_transcript_IDs_ecoli_UP_cutoff{cutoff}.txt", "w") as out_up, \
        open(f"Analysis_alle/transcripts_only/{cutoff}/alle_transcript_IDs_ecoli_DOWN_cutoff{cutoff}.txt", "w") as out_down, \
        open(f"Analysis_alle/transcripts_only/{cutoff}/alle_transcript_IDs_ecoli_all_cutoff{cutoff}.txt", "w") as out_all:
        for n in sorted(ecoli_SMESTs_UP):
            out_up.write(n + "\n")
            out_all.write(n + "\n")
        for n in sorted(ecoli_SMESTs_DOWN):
            out_down.write(n + "\n")
            out_all.write(n + "\n")
