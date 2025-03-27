import re
CELL_SIZE = 68
BOUNDARY = (round(-1 * CELL_SIZE/2,2),round(CELL_SIZE/2,2))
LEFT_BOUNDARY = BOUNDARY[0]
RIGHT_BOUNDARY = BOUNDARY[1]
NO_OF_BUCKET = 512
BUCKET_LENGTH = CELL_SIZE/ (NO_OF_BUCKET**(1/3))

def parsed_data(data):
        x_data = float(data[2])
        y_data = float(data[3])
        z_data = float(data[4])

        #x_cord
        # if x_data > CELL_SIZE/2:
        #     x_data = x_data - CELL_SIZE

        # #y_cord
        # if y_data > CELL_SIZE/2:
        #     y_data = y_data - CELL_SIZE

        # #z_cord
        # if z_data > CELL_SIZE/2:
        #     z_data = z_data - CELL_SIZE

        return {
            "atom_id" : data[0],
            "atom_class" : data[1],
            "atom_coordinate" : [x_data,y_data,z_data] 
        }

def get_points(path,num):
        file = open(path,"r")
        read_it = file.read()

        lines = read_it.splitlines()
        check_next = False
        initial = True
        flag = 0
        final_data = []
        one_snap = []

        for line in lines:
                if check_next:
                        NUMBER_OF_ATOMS = int(line.strip())
                        check_next = False
                        break

                if line == "ITEM: NUMBER OF ATOMS" and initial:
                        initial = False
                        check_next = True

        for i,line in enumerate(lines):
                line = line.strip()
                if "ITEM: ATOMS id type x y z" in line:
                        if flag:
                                final_data.append(one_snap)
                                one_snap = []
                        flag = 1
                match = re.match("([0-9]*)([ ]*)([0-9]*)([ ]*)([0-9]*\.*\d*)([ ]*)([0-9]*\.*\d*)([ ]*)([0-9]*\.*\d*)",line)
                if match:
                        if match.group(3) == num:

                                one = match.group(1)
                                two = match.group(3)
                                three = match.group(5)
                                four = match.group(7)
                                five = match.group(9)

                                val = parsed_data([one if one else 0,two if two else 0,three if three else 0,four if four else 0,five if five else 0])
                                one_snap.append(val)

        pts = [[a for a in sorted(individual_data,key=lambda x:int(x['atom_id']))] for individual_data in final_data]
        return pts

if __name__ == '__main__':
        import json
        with open("atom_info.json","r") as file:
            atoms = json.load(file) 
        for i,atom in enumerate(atoms):
                data = get_points('out.dump',str(i+1))
            
                out_file = open("{}.json".format(atom), "w")

                json.dump(data, out_file, indent = 6)
                  
                out_file.close()







