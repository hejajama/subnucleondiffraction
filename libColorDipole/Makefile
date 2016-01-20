#  All making is done in "obj" directory, just redirect execution there
.DEFAULT:
	$(MAKE) -C obj $(MAKEFLAGS) $(MAKECMDGOALS)


# Below is the list of all executable files which are produced, needed for distclean
EXEC_TARGETS=Test_Dipole

.PHONY: distclean
distclean:
		rm -f $(EXEC_TARGETS)
		rm -f libs/*
		$(MAKE) -C obj distclean
