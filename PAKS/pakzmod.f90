      MODULE DRAKCE8
        INTEGER*8 nwk8
        INTEGER*8, DIMENSION(:), ALLOCATABLE :: MAXA8
        INTEGER*4, DIMENSION(:), ALLOCATABLE :: IVRS
        INTEGER*4, DIMENSION(:), ALLOCATABLE :: ISK
        INTEGER*4, DIMENSION(:), ALLOCATABLE :: AIROWS
        INTEGER*4 TIPTACKANJA
      END MODULE
      
      MODULE FSIDENT
        INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: IDENT
      END MODULE
      
      MODULE PLAST3D
        REAL*8, DIMENSION(:), ALLOCATABLE :: PLAST
        REAL*8, DIMENSION(:), ALLOCATABLE :: PLAS1
        REAL*8, DIMENSION(:), ALLOCATABLE :: PLAS0
      END MODULE
      
      MODULE ELEMENTI
        INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: VTKELEMENTI
      END MODULE
!     ovo sluzi za citanje / pisanje mesto fajla ali samo za 1 slucaj koji je pravio problem      
      MODULE WSTAZK
        LOGICAL ALLOCIRANAMATRICA
        REAL*8, DIMENSION(:), ALLOCATABLE :: RTWRITE
      END MODULE
!     ovo je novi modul za cuvanje matrice K      
      MODULE MATRICA
        LOGICAL ALLOCIRANAMATRICAK
        REAL*8, DIMENSION(:), ALLOCATABLE :: ALSK
        REAL*8, DIMENSION(:), ALLOCATABLE :: ALSM
        REAL*8, DIMENSION(:), ALLOCATABLE :: ALSC
    END MODULE

      module STIFFNESS
        integer*8 :: stiff_n
        integer*8 :: nonzeros
        integer*8,dimension(:),allocatable :: rows
        integer*8,dimension(:),allocatable :: columns
        !ova druga 2 su privremena dok ne regulisem razliku izmedju integer 4 i 8
        integer*4,dimension(:),allocatable :: iirows
        integer*4,dimension(:),allocatable :: iicolumns
        integer*4,dimension(:),allocatable :: IMAXA
        integer*4,dimension(:),allocatable :: IMHT
      end module STIFFNESS
     