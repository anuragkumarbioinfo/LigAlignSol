
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
        self.oxygen = map(float,oxygen.split())
        self.pdb_coords = []
    def test(self, outfile):
        self.pdb_Read()
        self.calculate_COG()
        self.write_PDB(outfile)
    def pdb_Read(self):
        with open(self.pdbfile,"r") as reader:
            for lines in reader:
                if lines.startswith("ATOM"):
                    lines = map(float,lines.split()[5:8])
                    self.pdb_coords.append(lines)
    
    def calculate_COG(self):
        x1,y1,z1 = zip(*self.pdb_coords)
        self.pdb_COG = [sum(x1)/len(x1),sum(y1)/len(y1),sum(z1)/len(z1)]
        
        
    def write_PDB(self,outfile):
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
                        output.write(lines[:26]+'   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n'.format(new_x[counter], new_y[counter],new_z[counter], 0.00, 0.00))
                        counter += 1
        
