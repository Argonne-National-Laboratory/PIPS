lib_LIBRARIES = libhsl_ma86.a libfakemetis.a

modulesdir = modules
nodist_modules_HEADERS = \
	hsl_ma86_single.mod \
	hsl_ma86_double.mod \
	hsl_ma86_complex.mod \
	hsl_ma86_double_complex.mod

libhsl_ma86_a_SOURCES = \
	hsl_ma86s.f90 \
	hsl_ma86d.f90 \
	hsl_ma86c.f90 \
	hsl_ma86z.f90 \
	common.f \
	common90.f90 \
	sdeps.f \
	sdeps90.f90 \
	ddeps.f \
	ddeps90.f90 \
	cdeps.f \
	cdeps90.f90 \
	zdeps.f \
	zdeps90.f90 \
	hsl_ma86s_ciface.f90 \
	hsl_ma86d_ciface.f90 \
	hsl_ma86c_ciface.f90 \
	hsl_ma86z_ciface.f90 \
	hsl_mc68i_ciface.f90 \
	hsl_mc69s_ciface.f90 \
	hsl_mc69d_ciface.f90 \
	hsl_mc69c_ciface.f90 \
	hsl_mc69z_ciface.f90

libfakemetis_a_SOURCES = fakemetis.f

hsl_ma86s.$(OBJEXT): sdeps90.$(OBJEXT)
hsl_ma86d.$(OBJEXT): ddeps90.$(OBJEXT)
hsl_ma86c.$(OBJEXT): cdeps90.$(OBJEXT)
hsl_ma86z.$(OBJEXT): zdeps90.$(OBJEXT)
sdeps90.$(OBJEXT): common90.$(OBJEXT)
ddeps90.$(OBJEXT): common90.$(OBJEXT)
cdeps90.$(OBJEXT): common90.$(OBJEXT)
zdeps90.$(OBJEXT): common90.$(OBJEXT)
hsl_ma86s_ciface.$(OBJEXT): hsl_ma86s.$(OBJEXT)
hsl_ma86d_ciface.$(OBJEXT): hsl_ma86d.$(OBJEXT)
hsl_ma86c_ciface.$(OBJEXT): hsl_ma86c.$(OBJEXT)
hsl_ma86z_ciface.$(OBJEXT): hsl_ma86z.$(OBJEXT)
hsl_ma86i_ciface.$(OBJEXT): hsl_ma86i.$(OBJEXT)
hsl_ma86s_ciface.$(OBJEXT): hsl_ma86s.$(OBJEXT)
hsl_ma86d_ciface.$(OBJEXT): hsl_ma86d.$(OBJEXT)
hsl_ma86c_ciface.$(OBJEXT): hsl_ma86c.$(OBJEXT)
hsl_ma86z_ciface.$(OBJEXT): hsl_ma86z.$(OBJEXT)
