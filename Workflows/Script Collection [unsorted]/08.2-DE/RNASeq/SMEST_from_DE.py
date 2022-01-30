import os
import json

cutoff = 1.9 

output_dir = os.getcwd()

DE_Analysis_dir = "/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08.2-DE/RNASeq/01-DE/Analysis/"
polyIC_spec_file = "de_polyIC_specific_transcriptsonly.json"
ecoli_spec_file= "de_ecoli_specific_transcriptsonly.json"
polyIC_all_file="de_unique_transcriptsonly2plus_polyic.json"
ecoli_all_file="de_unique_transcriptsonly2plus_ecoli.json"

polyIC_spec_dict, ecoli_spec_dict = json.load(open(DE_Analysis_dir+polyIC_spec_file)), json.load(open(DE_Analysis_dir+ecoli_spec_file))
polyIC_spec_up_dict, ecoli_spec_up_dict = polyIC_spec_dict.get("UP"), ecoli_spec_dict.get("UP")
polyIC_spec_down_dict, ecoli_spec_down_dict = polyIC_spec_dict.get("DOWN"), ecoli_spec_dict.get("DOWN")

polyic_spec_down, polyic_spec_up, ecoli_spec_down, ecoli_spec_up, \
polyic_all_down, polyic_all_up, ecoli_all_down, ecoli_all_up = [list() for x in range(0, 8)]

for range, transcripts in polyIC_spec_down_dict.items():
    for transcript, information in transcripts.items():
        if abs(information.get("log2FoldChange")) < cutoff:
                continue
        polyic_spec_down.append(transcript)
for range, transcripts in polyIC_spec_up_dict.items():
    for transcript, information in transcripts.items():
        if abs(information.get("log2FoldChange")) < cutoff:
                continue
        polyic_spec_up.append(transcript)
for range, transcripts in ecoli_spec_down_dict.items():
    for transcript, information in transcripts.items():
        if abs(information.get("log2FoldChange")) < cutoff:
                continue
        ecoli_spec_down.append(transcript)
for range, transcripts in ecoli_spec_up_dict.items():
    for transcript, information in transcripts.items():
        if abs(information.get("log2FoldChange")) < cutoff:
                continue
        ecoli_spec_up.append(transcript)

polyIC_all_dict, ecoli_all_dict = json.load(open(DE_Analysis_dir+polyIC_all_file)), json.load(open(DE_Analysis_dir+ecoli_all_file))
polyIC_all_up_dict, ecoli_all_up_dict = polyIC_all_dict.get("UP"), ecoli_all_dict.get("UP")
polyIC_all_down_dict, ecoli_all_down_dict = polyIC_all_dict.get("DOWN"), ecoli_all_dict.get("DOWN")

for range, transcripts in polyIC_all_down_dict.items():
    for transcript in transcripts:
        if not transcript.endswith("1.1"):
                continue
        polyic_all_down.append(transcript)
for range, transcripts in polyIC_all_up_dict.items():
    for transcript in transcripts:
        if not transcript.endswith("1.1"):
                continue
        polyic_all_up.append(transcript)
for range, transcripts in ecoli_all_down_dict.items():
    for transcript in transcripts:
        if not transcript.endswith("1.1"):
                continue
        ecoli_all_down.append(transcript)
for range, transcripts in ecoli_all_up_dict.items():
    for transcript in transcripts:
        if not transcript.endswith("1.1"):
                continue
        ecoli_all_up.append(transcript)

os.chdir("DE_unique_Transcripts")
for typus in ["ecoli", "polyic"]:
    try:
        os.mkdir(typus)
    except:
        pass
    for updown in ["up", "down"]:
        try:
            os.mkdir(f"{typus}/{updown}")
        except:
            pass
        for allspec in ["all", "spec"]:
            with open(f"{typus}/{updown}/{typus}_{allspec}_{updown}.txt", "w") as out:
                for smest in locals()[f"{typus}_{allspec}_{updown}"]:
                    out.write(f"{smest}\n")
            with open(f"{typus}_{allspec}_{updown}.txt", "w") as out:
                for smest in locals()[f"{typus}_{allspec}_{updown}"]:
                    out.write(f"{smest}\n")

for typus in ["ecoli", "polyic"]:
        for allspec in ["all", "spec"]:
                with open(f"{typus}/{typus}_{allspec}_upanddown.txt", "w") as out:
                            for smest in locals()[f"{typus}_{allspec}_up"]+locals()[f"{typus}_{allspec}_down"]:
                                out.write(f"{smest}\n")
                with open(f"{typus}_{allspec}_upanddown.txt", "w") as out:
                            for smest in locals()[f"{typus}_{allspec}_up"]+locals()[f"{typus}_{allspec}_down"]:
                                out.write(f"{smest}\n")
