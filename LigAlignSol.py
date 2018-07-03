
class Geometry:
    """
    LigAlignSol.py is a python script which superimpose or align Ligand coordinates to the oxygen atom Water coordinates(or SOL as per Gromacs naming convention) or heavy atom of any other coordinate system.  
    This script may help to place the ligand/substrate at the position of oxygen atom of water(SOL) or heavy atom of any coordinate.
    
    LigAlignSol v0.1
    
    Input: 1. Ligand pdb file required.
           2. Oxygen atom coordinate(x,y,z) required.
    Output: Ligand pdb superimposed with oxygen atom
    
    Method: 1. Read ligand pdb file
            2. Center Of Geometry(C.O.G) of ligand pdb file calculated
            3. Write superimposed ligand output pdb file as output.pdb 
    """
    def __init__(self,pdbfile,oxygen):
        self.pdbfile = pdbfile
        if type(oxygen) == str: oxygen = oxygen.split()
        self.oxygen = map(float,oxygen)
        self.pdb_coords = []
        self.distance = lambda a,b: sum([(a[i] - b[i])**2 for i in range(len(a))]) ** 0.5
        self.savedistance = {}
    def runscript(self, outfile):
        self.pdb_Read()
        self.calculate_COG()
        self.write_PDB(outfile)
        
    def pdb_Read(self):
        with open(self.pdbfile,"r") as reader:
            for lines in reader:
                if lines.startswith("ATOM"):
                    lines = map(float,lines.split()[5:8])
                    self.pdb_coords.append(lines)
    
    def calculate_COG(self, final=False):
        #print "[COG]: ", final
        if not final: 
            x1,y1,z1 = zip(*self.pdb_coords)
            self.pdb_COG = [sum(x1)/len(x1),sum(y1)/len(y1),sum(z1)/len(z1)]
        else: 
            x1,y1,z1 = zip(*final)
            return([sum(x1)/len(x1),sum(y1)/len(y1),sum(z1)/len(z1)])
    def write_PDB(self,outfile):
        self.otf = []
        Nx, Ny, Nz = self.oxygen
        Ox, Oy, Oz = self.pdb_COG
        new_x = []
        new_y = []
        new_z = []
        def XYZN(I,N,O):
            return (I + N)- O
        for i,j,k in self.pdb_coords:
            new_x.append(XYZN(i,Nx, Ox))
            new_y.append(XYZN(j,Ny, Oy))
            new_z.append(XYZN(k,Nz, Oz))
        counter = 0
        with open(self.pdbfile,"r") as p1_g:
            with open(outfile,"w") as output:
                for lines in p1_g:
                    if lines.startswith("ATOM"):
                        output.write(lines[:26]+'   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n'.format(new_x[counter]+11.00, new_y[counter],new_z[counter], 0.00, 0.00))
                        self.otf.append([new_x[counter]+11.00, new_y[counter],new_z[counter]])
                        counter += 1
    
    def nearest_atom(self, coordinate,protein):
        self.coordinate = map(float,coordinate)
        self.coords = {}
        distance_list = []
        with open(protein) as pp:
            for lines in pp:
                if lines.startswith("ATOM"):
                    linesZ = map(float,lines.split()[5:8])
                    name = "_".join(lines.split()[1:5])
                    self.coords[name] = linesZ
                    self.savedistance[name] = self.distance(linesZ,self.coordinate)
        rr = min(self.savedistance, key=self.savedistance.get)
        #print self.calculate_COG(self.otf)
        print "\n\nDistance between {} and ligand is = {:5.5} Angstrom".format(rr, self.distance(self.calculate_COG(self.otf),self.coords[rr] ))
        return rr, self.savedistance[rr]



import sys
#python LigAlignSol.py ACX.pdb tail 94.935  50.098  25.530 head 74.398  49.079  36.679 protein.pdb output.pdb
a = sys.argv
pdb = a[1]
protein = a[-2]
xyz = a[2:5]
pxyz = a[5:8]
output = a[-1]
query = Geometry(pdb, xyz)
query.runscript(output)
x,d = query.nearest_atom(pxyz,protein)
print "Nearest Neighbour found: {} at {:5.5} Angstrom\n\n".format(x,d)
#print query.distance([1,1,1],[0,0,0])
