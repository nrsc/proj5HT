PROJECT_HOME := $$HOME/proj/proj5HT
PROJECT_DEN :=	$$HOME/proj/proj5HT/den/

NHP_DATA := /allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/NHP_expts
HUMAN_DATA := /allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/Human_expts
Seq_DATA := /allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Macaque_NCBI
sert_repo := /allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/manuscript_code/serotonin
ROOKERY := /allen/programs/celltypes/workgroups/hct/SawchukS/rookery

SYMLINKS := \
	$(PROJECT_HOME):$(ROOKERY) \
	$(PROJECT_DEN):$(NHP_DATA) \
	$(PROJECT_DEN):$(HUMAN_DATA) \
	$(PROJECT_DEN):$(Seq_DATA) \
	$(PROJECT_DEN):$(sert_repo)

.PHONY: all symlinks
all: symlinks


symlinks:
	@set -e; \
	OSNAME=$$(uname 2>/dev/null || echo Windows); \
	for pair in $(SYMLINKS); do \
		dest=$${pair%%:*}; \
		src=$${pair##*:}; \
		linkname="$${dest}/$$(basename $$src)"; \
		if [ -L "$$linkname" ]; then \
			echo "Symlink already exists: $$linkname"; \
		elif [ -e "$$linkname" ]; then \
			echo "Error: $$linkname exists but is not a symlink"; \
			exit 1; \
		else \
			echo "Creating symlink: $$linkname -> $$src"; \
			if [ "$$OSNAME" = "Linux" ] || [ "$$OSNAME" = "Darwin" ]; then \
				ln -s "$$src" "$$linkname"; \
			elif [ "$$OSNAME" = "Windows" ]; then \
				powershell -Command "New-Item -ItemType SymbolicLink -Path '$$linkname' -Target '$$src'" || exit 1; \
			else \
				echo "Unsupported OS: $$OSNAME"; \
				exit 1; \
			fi; \
		fi; \
	done





