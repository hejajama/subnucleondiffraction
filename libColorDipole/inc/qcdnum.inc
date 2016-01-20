
      parameter ( mxg0   = 5       )   !max # x-grids
      parameter ( mxx0   = 1000     )   !max # x-gridpoints
      parameter ( mqq0   = 1000     )   !max # Q2-gridpoints
      parameter ( mpt0   = 5000    )   !max # interpolations in one call
      parameter ( miw0   = 20      )   !max # info words in weight store
      parameter ( mbf0   = 20      )   !max # scratch buffers fast engine
      parameter ( nwf0   = 2500000 )   !tot # of words in pdf store

C--   --------- do not touch the values below ---------------------------

      parameter ( iom0   = 3      )    !max spline order in x or Q2
      parameter ( mnf0   = 12     )    !# q + qbar flavors
      parameter ( mbp0   = 3      )    !# of status words
      parameter ( nwrk0  = 4      )    !# working tables (scratch space)
      parameter ( mpp0   = 7      )    !# Pij tables (QQ,QG,GQ,GG,PL,MI,VA)
      parameter ( maa0   = 5      )    !# Aij tables
      parameter ( mset0  = 9      )    !# pdf sets (unpol,pol,frag,custom,external)
