from phenix.programs import homology
from mmtbx.nci import hbond
import time
import mmtbx
import mmtbx.alignment
from libtbx import easy_run
from libtbx import group_args
import iotbx.pdb
import iotbx.pdb.fetch
import iotbx.phil
from libtbx import easy_pickle
from libtbx.utils import Sorry
from iotbx import bioinformatics
import iotbx.bioinformatics.pdb_info
from iotbx.bioinformatics import local_blast
from iotbx.bioinformatics.structure import summarize_blast_output
from cctbx.array_family import flex


def prepare_hydrogen_restraints(pdb_file):
  params = homology.get_default_params()
  params.num_of_best_pdb = 1
  pdb_code_ref = pdb_file[0:4]
  res = homology.file_perfect_pair(pdb_file, params)
  for r in res :
    data_type_A = False
    print (r.chain_ref,r.match[0].chain_id,r.match[0])
    pdb_code_ref = pdb_file[0:4]
    pdb_code = r.match[0].pdb_code
    chain_id = r.match[0].chain_id
    # easy_run.call("phenix.fetch_pdb {0}".format(pdb_code))
    pdb_inp = iotbx.pdb.input(pdb_code+".pdb")
    data_type_A =pdb_inp.get_experiment_type()
    if data_type_A == False:continue
    he = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
    hx = iotbx.pdb.input(file_name=pdb_code + ".pdb").construct_hierarchy()
    sel_e = he.atom_selection_cache().selection("protein and chain %s"%(r.chain_ref))
    he = he.select(sel_e)
    sel_x = hx.atom_selection_cache().selection("protein and chain %s" % (chain_id))
    hx = hx.select(sel_x)
    chain_A = (he.only_chain())
    chain_B = (hx.only_chain())
    #phil_obj = align_tow_chain(chain_A,chain_B,phil_obj)
def align_tow_chain(chain_A,chain_B,phil_obj):
  se = chain_A.as_padded_sequence()
  sx = chain_B.as_padded_sequence()
  align = mmtbx.alignment.align(
    seq_a = se,
    seq_b = sx)
  alignment = align.extract_alignment()
  matches = alignment.identity_matches()
  a_seq = alignment.a
  b_seq = alignment.b
  ai_seq = alignment.i_seqs_a
  bi_seq = alignment.i_seqs_b
  print (dir(alignment))
  print (a_seq)
  print (matches)
  print (b_seq)
  print (ai_seq)
  print (bi_seq)
  return phil_obj
if __name__ == '__main__':
    start = time.time()
    prepare_hydrogen_restraints('6d9h.pdb')
    end = time.time()
    time_cost = (end - start)
    print "it cost % seconds" % time_cost
