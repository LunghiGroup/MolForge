.SUFFIXES:

targets = CF abc2cell cell2abc Bispectrums C62LJ Distribution Distribution2D FindMols Correlation Eckart getFC Mat2Euler get_topo FT1d SuperCellGen fitcharges 

all : $(addprefix Tool_,$(targets)) 

Tool_%: 
	cd $(TOOLS_DIR)/$(subst Tool_,,$@) && $(MAKE) $(subst Tool_,,$@)

clean : $(addprefix Clean_,$(targets)) 

Clean_%: 
	cd $(TOOLS_DIR)/$(subst Clean_,,$@) && $(MAKE) clean
