# TODO::
## The implementation of everything in this file is hacky beyond belief


from Zernike3 import *
from Bio.PDB import PDBList
from os import path, close
import tempfile
import shutil

class Temp_PDB_Downloader:
    pdbl = PDBList()
    pdb = None
    temp_dir = None
    temp_file = None
    def download_pdb(self, pdb):
        self.temp_dir = tempfile.mkdtemp()
        # This does not return the filename it returns something weird
        self.temp_file = pdbl.retrieve_pdb_file(pdb, pdir = self.temp_dir, file_format = '.pdb') 
        return(filename)
    def remove_pdb(self):
        close(temp_file)
        rmtree(self.temp_dir)
        return(pdb_file_path)


def voxelize_pdb(pdb, resolution, maxext=0.6, splatfrac=0.28, save_dir=None, grid_dir=None):
    tpd = None
    if not pdb.endswith('.pdb'):
        pdb_id = pdb
        tpd = Temp_PDB_Downloader()
        pdb = tpd.download_pdb(pdb)
    else:
        pdb_id = path.splitext(path.basename(pdb))[0]
    vox = Voxels()
    vox.SetResolution(resolution)
#     vox.SetMaxExtent(maxext)
#     vox.SetSplatFrac(splatfrac)
    vox.LoadStructure(pdb)
    if not save_dir is None:
        vox.SaveVoxels('{0}/{1}.pkl'.format(path.normpath(save_dir), pdb_id))
    if not grid_dir is None:
        vox.Grid2DX('{0}/{1}.dx'.format(path.normpath(grid_dir), pdb_id)) 
    if tpd is not None:
        tpd.remove_pdb()
    return(vox.voxels)

def voxelize_pdb_list(pdb_path_list, resolution, maxext=0.6, splatfrac=0.28, pickle_file=None):
    out = dict()
    for pdb_path in pdb_path_list:
        fname = path.splitext(path.basename(pdb_path))[0]
        pdbvox = voxelize_pdb(pdb_path, resolution, maxext, splatfrac)
        out.update({fname : pdbvox})
    if not pickle_file is None:
        fname = open(pickle_file, 'wb')
        pickle.dump(fname, out)
    return(out)
                              
                              
if __name__ == '__main__':
    res = voxelize_pdb('ZTEST/1gww.pdb', resolution=1000)
#     pdb_path_list = ['ZTEST/1gww.pdb', 'ZTEST/1ej1.pdb']
#     res = voxelize_pdb_list(pdb_path_list, resolution=2)
    print(res)