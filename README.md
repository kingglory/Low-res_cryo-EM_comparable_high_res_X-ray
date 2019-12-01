# Low-res_cryo-EM_comparable_high_res_X-ray
Implementation stages:

1. Write a function that takes input atomic model (pdb.hierarchy object or mmtbx.model object, TBD) and returns
geometry.edits (as Phil object) that defines H-bonds. To do this, for each chain of low-res model find matching chain in
high-res model. For this use underlying code of phenix.homology. Once each chain in low-res model is associated with
the corresponding chain from high-res model (or none, if no match found), the function:
  a. Finds all H-bonds in high-res model.
  b. Identifies corresponding (matching) atoms in low-res model. For this one may need to align both chains, then for
     each set of atoms making H bond in high-res model find corresponding set of atoms in the low-res model. Write this
     as a separate function that takes two chains (chain object of PDB hierarchy) and returns H bond definitions as
     geometry.edits Phil object.
  c. Dump geometry.edits Phil object into a Phil parameter file.
  
2. Once the function described in #2 is available:
  a. Find out how many cryo-EM models have matching higher-resolution X-ray models.
     i. For each cryo-EM model find out percentage of its atoms that have a match in a high-res X-ray model.
     ii. Come up with a way to present these statistics efficiently (a graph, histogram, scatter plot, etc).
  b. Using results from “a”, select 5-10 cryo-EM models for which high-resolution matching models are available.
  c. These cryo-EM models must possess poor validation metrics and poor H bonding parameters.
  d. Each chain of cryo-EM model must have fully matching chain from a high-res model.
  e. High-res model must have all acceptable validation metrics. For this consider exploring several choices of matching
     high-res models that phenix.homology can produce.
  f. Once “cryo-EM <> high-res” pairs are selected:
     i. Generate H bond restraints (geometry.edits Phil file) using the tool from Stage 1.
     ii. Run real-space refinement of selected cryo-EM models using three different modes:
       • Refine with all default settings
       • Refine using high-res model as a reference model (XXX)
       • Refine using H-bond restraints derived from high-res model
     iii. Analyze quality of refined models (XXX).
     
     
     For stage 1 function: hydrogen_bond_restraints.py
     For stage 2 :6h3n.pdb 100% matching high resolution
                  5m6s.pdb 100%
                  6n2p.pdb 100%
                  6h3l.pdb 100%
                  601m.pdb 97.86%
                  601l.pdb 97.64%
                  6r93.pdb 95.08%
