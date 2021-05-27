!
!
!****************************************************************************
!
!  SUBROUTINE: LecSampler
!
!  PURPOSE:  Read input sampler data.
!
!****************************************************************************
subroutine LecSampler()
!
!
!
!	Variables globales
	use fslamGlobals_structures
	use fslamGlobals_shared
	use fslamGlobals_private
!
!
	implicit double precision (a-h,o-z)
!
!
!
!	Abrimos ficheros de entrada de datos de control
	open(100,file=fname_sampler_dat,status='old',form='formatted')
!
	write(6,'("Reading sampling file...",/)')
!
!	Entrada de parametros de input.dat:
    read(100,*) rMinAR      !Minimum Antecedent Rainfall
!
    read(100,*) rMaxAR      !Maximum Antecedent Rainfall
!
    read(100,*) iNumAR      !Number of Atecedent Rainfall values
!
    read(100,*) rMinER      !Minimum Event Rainfall
!
    read(100,*) rMaxER      !Maximum Event Rainfall
!
    read(100,*) iNumER      !Number of Event Rainfall values
!
    read(100,*) rThres      !PoF threshold
!
!	Cerramos el fichero
    close(100)
!
!	Dynamic array
	ALLOCATE(samplingRes(iNumAR + 1,iNumER + 1))			        !Zones delimitation grid

!
!    
!
end subroutine
