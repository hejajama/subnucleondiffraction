
      logical mark9
      common /fmark9/ mark9(0:mxx0,0:mqq0+7)
      common /fmesh9/ iy19(mpt0),iy29(mpt0),iz19(mpt0),iz29(mpt0)
      common /flist9/ xlst9(mpt0),qlst9(mpt0),
     +                ylst9(mpt0),tlst9(mpt0),
     +                ixqfyt9(mpt0),nxq9,nyt9
      common /sparse/ nzlist9,izlist9(mqq0),nyslist9(mqq0),
     +                iyslist9(mxx0,mqq0),nydlist9(mqq0),
     +                iydlist9(mxx0,mqq0),isparse9(mbf0) 

