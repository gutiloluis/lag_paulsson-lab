import csv
import numpy as np
import re

def parse_bit_command(bit_command):
    # Use regular expressions to extract the integer from the command
    match = re.match(r'BIT(\d+)(ON|OFF)', bit_command)
    if match:
        # Extract the integer and the state (ON or OFF)
        bit_number = int(match.group(1))
        state = 1 if match.group(2) == 'ON' else 0
        return bit_number, state
    else:
        # Return None if the command doesn't match the expected format
        return None
        
def convert_to_binary(input_string):
    # Split the input string by '>' and '<' to extract individual commands
    commands = re.split('[><]', input_string)
    if commands[1] == 'UNS3':
        sign = '+'
    else:
        sign = '-'
    commands = [cmd for cmd in commands if 'BIT' in cmd]

    # Initialize a list to store binary values
    
    # Iterate through each command and determine if it turns a bit on or off
    if len(commands) == 30:
        bit_list = [None]*30
        for cmd in commands:
            bit_number, bit_state = parse_bit_command(cmd)
            bit_list[bit_number] = bit_state
            binary_string = ''.join(map(str, bit_list))
    else:
        binary_string = ''
    return sign, binary_string

cigar_dict = {}
# with open(snakemake.input[0], "r") as infile:
#     for line in infile:
#         data = line.split("\t")
#         read_id = data[0].split(" ")[0]
#         cigar_dict[read_id] = (
#             int(data[7]),
#             int(data[8]),
#             data[5],
#             data[15].split(":")[-1][:-1],
#             )
with open(snakemake.input[0], "r") as infile:
    for line in infile:
        data = line.split("\t")
        read_id = data[0]
        cigar_dict[read_id] = (
            int(data[5]),
            int(data[6]),
            data[8],
            data[19].split(":")[-1][:-1],
            )        
    barcode_dict = {}
    for key in cigar_dict.keys():
        cigar = cigar_dict[key]
        if ("UNS3" in cigar[2]) and ("UNS9" in cigar[2]):

            # print("valid read")
            sign, barcode = convert_to_binary(cigar[2])
            barcode_dict[key] = barcode
            cigar_dict[key] = tuple(sign)+cigar_dict[key]
    barcode_dict = [
    {"readname": key, "barcode": val, "cigar": cigar_dict[key]}
    for key, val in barcode_dict.items()]


keys = ["readname", "barcode", "cigar"]
with open(snakemake.output[0], "w") as outfile:
    dict_writer = csv.DictWriter(outfile, keys, delimiter="\t")
    dict_writer.writeheader()
    dict_writer.writerows(barcode_dict)





# import csv
# import numpy as np

# cigar_dict = {}
# with open(snakemake.input[0], "r") as infile:

#     for line in infile:
#         data = line.split("\t")
#         read_id = data[0].split(" ")[0]
#         if ">" in data[5]:
#             cigar_dict[read_id] = (
#                 "+",
#                 int(data[7]),
#                 int(data[8]),
#                 data[5],
#                 data[15].split(":")[-1][:-1],
#             )
#         else:
#             cigar_dict[read_id] = (
#                 "-",
#                 int(data[7]),
#                 int(data[8]),
#                 data[5],
#                 data[15].split(":")[-1][:-1],
#             )

# barcode_dict = {}
# for key in cigar_dict.keys():
#     cigar = cigar_dict[key]
#     if "Variant" in cigar[3] and "SPACER2" in cigar[3]:
#         if cigar[0] == "-":
#             barcode = cigar[3].split("<")
#             barcode = barcode[::-1]
#         elif cigar[0] == "+":
#             barcode = cigar[3].split(">")
#         barcode = barcode[:-1]
#         barcode = (
#             np.array(["ON" in item for item in barcode if "BIT" in item])
#             .astype(int)
#             .astype(str)
#             .tolist()
#         )
#         barcode = "".join(barcode)
#         barcode_dict[key] = barcode

# barcode_dict = [
#     {"readname": key, "barcode": val, "cigar": cigar_dict[key]}
#     for key, val in barcode_dict.items()
# ]

# keys = ["readname", "barcode", "cigar"]
# with open(snakemake.output[0], "w") as outfile:
#     dict_writer = csv.DictWriter(outfile, keys, delimiter="\t")
#     dict_writer.writeheader()
#     dict_writer.writerows(barcode_dict)
