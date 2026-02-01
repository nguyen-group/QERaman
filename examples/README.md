- Example 1: The Raman is calculated for the MoS2 monolayer, in which we set a dense k-mesh of 48x48x1 in the SCF file.
- Example 2: The Raman is calculated for the MoS2 monolayer, in which we set a coarse k-mesh of 12x12x1 in the SCF file, and a dense k-mesh of 48x48x1 in the NSCF file

Note: Example 2 shows a significant reduction in time (7.5-fold) compared with example 1. The workflow is: SCF (coarse k-mesh) --> PH (only fildvscf = 'dvscf') --> NSCF (dense k-mesh) --> BANDS_MAT (dipol vector for dense k-mesh) --> PH_MAT (el-ph for a dense k-mesh by using previously saved DeltaVscf in fildvscf).
