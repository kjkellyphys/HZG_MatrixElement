      double complex function A1fggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (100)
      include 'constants.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      double complex BSYA1fggpppp
      integer e1,p2,p3,e4
 
      A1fggpppp=BSYA1fggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
