#---------------------------------------------------------------------------------------------
# Parse the cabling map obtained with DTCCablingMapProducer_dump.py, 
# and sort it the same way as expected by the unpacker code.
# No info on module type, this is just to compare with info from Phase2TrackerCabling. 
#---------------------------------------------------------------------------------------------

import os, pdb

def process_file_main(input_file_main, output_file):

    # Read the main input file and collect the rows
    rows = []
    with open(input_file_main, 'r') as infile:
        # Skip first N lines 
        for i in range(16):
            next(infile); 
        
        for line in infile:
            # skip the last part of the file  
            if '%MSG' in line: 
                break
            # Split the line by commas to get the three integers
            numbers = line.split(',')
            if len(numbers) > 1: # Skip empty lines
                det_id = int(numbers[0])
                dtc_id = int(numbers[1])
                dtc_ch = int(numbers[2])

                # Append a tuple of the four integers to the rows list
                rows.append((det_id, dtc_id, dtc_ch))
    
    # Sort the rows based on the third integer, then by the second integer if third integers are identical
    rows.sort(key=lambda x: (x[1], x[2]))
    
    # Write the sorted rows to the output file
    with open(output_filename, 'w') as outfile:
        for aa, bb, cc  in rows:
            
            # Write the result to the output file
            outfile.write(f'{aa}  {bb}  {cc}  \n')

#---------- Main ----------

input_filename_main = 'dump_TrackerDetToDTCELinkCablingMap_T33_v2.txt' 
output_filename = './tkLayoutCabling_DetToDTCElinkMap_parsed.txt'

if os.path.exists(output_filename):
  os.remove(output_filename)

process_file_main(input_filename_main,
                  output_filename)

print("File created:", output_filename)
