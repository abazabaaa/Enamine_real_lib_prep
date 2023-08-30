#PARENT_WORKING_DIRECTORY = '/cbica/home/grahamth/JOBS'
#PARENT_JOB_NAME = 'dopamine_dude_actives'
#PARENT_MOLECULE_DB = '/cbica/home/grahamth/dopamine_test1/actives_final.sdf'
#PARENT_MOLECULE_DB_EXT = '.sdf'
#PARENT_OUTPUT_DB_EXT = '.sdf'


# n_jobs = 20
# smiles_col = 'Smile'
# catalog_id_col = 'CatalogID'
# canonical_id_col = 'ID_Index'
# dataset_format = 'parquet'

RUN_NAME = 'decoys'

# SGE_CPU_PARAM = 'qsub -l h_vmem=8G'

SGE_CPU_PARAM = 'sbatch'

SGE_GPU_PARAM = 'qsub -l gpu -l h_vmem=24G'

# VENDOR_MOLECULE_DIR = '/cbica/home/grahamth/vendor_databases/enamine'
#job1
JOB_NAME = 'actives'
WORKING_DIR = '/cbica/home/grahamth/sweetnothings/various_tests'
MOLECULE_DB = '/cbica/home/grahamth/dopamine_test1/decoys_final.sdf'
MOLECULE_DB_EXT = '.sdf'
OUTPUT_DB_EXT = '.pdb'

#subjob1
SUBJOB1_CPU_CORES = 5
OBABEL_SDF_TO_PDB = 'obabel {input_file_path} -O {output_file_path}'
SUBJOB1_EXE_CMD = OBABEL_SDF_TO_PDB
INPUT_DB_EXT = '.sdf'
SUBJOB1INPUT_DB_EXT = '.sdf'
SUBJOB1OUTPUT_DB_EXT = '.pdb'
OBABEL_DEPENDENCY = '''\
export PATH=$PATH:/cbica/home/grahamth/openbabel3/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cbica/home/grahamth/openbabel3/lib
'''
SUBJOB1_DEPENDENCY = OBABEL_DEPENDENCY

SUBJOB2_NAME = 'convertpdb'
SUBJOB2_CPU_CORES = 5
SUBJOB2INPUT_DB_EXT = '.pdb'
SUBJOB2OUTPUT_DB_EXT = '.pdbqt'
PREPARE_PDBQTLIGAND = '/cbica/home/grahamth/autodock/ADFRsuite_x86_64Linux_1.0/buildtest/bin/prepare_ligand -v -l {input_file_path} -o {output_file_path}'
SUBJOB2_EXE_CMD = PREPARE_PDBQTLIGAND
SUBJOB2_DEPENDENCY = 'pwd'
#subjob2
SUBJOB_NAME2 = 'dock_pdbqt'

SUBJOB3_GPU = 10
RECEPTOR_PATH = '/cbica/home/grahamth/autodock/dopamine/dude_d3/pocket2_fixer_moreatoms/rigidReceptor.maps.fld'
AUTODOCK_GPU_128WI = '''\
/cbica/home/grahamth/autodock/AutoDock-GPU/bin/autodock_gpu_128wi \
		-ffile /cbica/home/grahamth/autodock/dopamine/dude_d3/pocket2_fixer_moreatoms/rigidReceptor.maps.fld \
		-lfile {input_file_path} \
		-resnam {output_file_path} \
		-nrun 10
'''
SUBJOB3_EXE_CMD = AUTODOCK_GPU_128WI
SUBJOB3_DEPENDENCY = 'nvidia-smi'
SUBJOB3INPUT_DB_EXT = '.pdbqt'
SUBJOB3OUTPUT_DB_EXT = '.top_poses'

#CPU_CORES = 500
#LIGANDS_PER_GPU = 1000

#OPEN_EYE_LICENSE = '/cbica/home/grahamth/oe_license.txt'






#PREPARE_PDBQTLIGAND = templates.PREPARE_LIGAND_SH



#OMEGA_CMD = '/cbica/home/grahamth/openeyeapp_2020_2/openeye/arch/redhat-RHEL7-x64/omega/4.1.0.0/oeomega pose -in {input_file_path} -out {output_file_path} -maxconfs 800'
#FRED_CMD = '/cbica/home/grahamth/openeyeapp_2020_2/openeye/arch/redhat-RHEL7-x64/oedocking/fred -receptor /cbica/home/grahamth/test_code/dopamine/pocket2_nowater.oedu -in {input_file_path} -out {output_file_path} -hitlist_size 0'
#OPTLIG_CMD = '/cbica/home/grahamth/openeyeapp_2020_2/openeye/arch/redhat-RHEL7-x64/szybki/2.3.0.0/optligandindu -du /cbica/home/grahamth/test_code/dopamine/pocket2_nowater.oedu -in {input_file_path} -out {output_file_path}'
#SZYBKI_SCORE_CMD = '/cbica/home/grahamth/openeyeapp_2020_2/openeye/arch/redhat-RHEL7-x64/szybki/2.3.0.0/szybki -protein /cbica/home/grahamth/test_code/dopamine/_AC__DU__.pdb -in {input_file_path} -out {output_file_path} -sdtag all -optGeometry sp -protein_elec PB -strip_water'
#HYBRID_CMD = '/cbica/home/grahamth/openeyeapp_2020_2/openeye/arch/redhat-RHEL7-x64/oedocking/hybrid -receptor /cbica/home/grahamth/test_code/dopamine/pocket2_nowater.oedu -in {input_file_path} -out {output_file_path} -hitlist_size 0'
#MOLCHARGE_CMD = '/cbica/home/grahamth/openeyeapp_2020_2/openeye/bin/molcharge -method am1bccsym -in {input_file_path} -out {output_file_path}'
#OBABEL_SDF_TO_MULTIPDB = 'obabel {input_file_path} -O {output_file_path} -m'




