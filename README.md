# SARS-COV-nanobody-screening
Screening of Nanobodies against SARS-CoV-2 and Identification of Active Residues
# Description
Nanobodies, also known as VHHs or Nbs, are single-domain antibodies derived from camelid heavy-chain antibodies. Nanobodies have many unique properties that make them good candidates for the treatment of viruses. However, there is a lack of in silico screening techniques to promptly identify a large number of nanobodies. In  this project, three-step screening strategy iteratively in Silico was applied to acquire effective anti-SARS-CoV-2 nanobodies. In the strategy, we use NeighborSearch to identify the active residues of nanobodies and antigen. The output file listed the IDs of the active residues.
# Usage
  ## Download data
   The required data can be downloaded from the respective database:
   
       * '8GZ5' is the PDBID，it can be replaced with other anti-SARS-Cov-2 nanobodies.
       * Please note that the ID needs to be capitalized when downloading using the following link.
       · Download the sample PDB file: http://www.nanolas.cloud/download/pdb/8GZ5.pdb   
  
  ## Identify Active Residues
  ### Active Residues of nanobodies
   Download the code: 
   
      · Download Identify_Active_Residues/toparquest.py
      
   Make sure to put **toparquest.py** and **8gz5.pdb** in the same directory.
   ### Active Residues of antigens
   Download the code: 
   
     · Download Identify_Active_Residues/GetClusproe.py

  Download the sample PDB:
  
     * The modl.pdb is a conformation of the Spike protein (antigen) bound with a nanobody.
     * It is used to identify the active residues of the Spike protein. The file does not contain COMPND sections.
     · Download Sample PDB/modl.pdb

  Make sure to put **GetClusproe.py** and **modl.pdb** in the same directory.
    
