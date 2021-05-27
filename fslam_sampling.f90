!  fslam_sampling.f90 
!
!  FUNCTIONS:
!	fslam_sampling      - Entry point of console application.
!
!
!	Read dem
!	Read data input
!	Read soil zones raster
!	Read soil properties
!	Read landuse raster
!	Read landuse properties
!	Fill sinks
!	Compute slopes
!	Compute flow accumulation
!   Compute stability
!

!****************************************************************************
!
!  PROGRAM: fslam_sampling
!
!  PURPOSE:  Fast Shallow Landslide Assessment Model
!
!****************************************************************************

program fslam_sampling
!
!	Variables globales
	use fslamGlobals_structures
	use fslamGlobals_shared
	use fslamGlobals_private
!
	implicit double precision (a-h,o-z)
!
!   Local variables
    integer :: i
    integer :: j
    real*8  :: numNullCells
    real*8  :: numNotNullCells
!
!
!	Printout
	write(6,'("*******************************************************")')
	write(6,'("Probabilistic Safety Factor model for shallow landslide")')
   	write(6,'("                    SAMPLING VERSION                   ")')
	write(6,'("*******************************************************",/)')
!
!   Read input filenames
    call LecFilenames()
!
!	Log file
	open(unit=100,file=(trim(fname_res) // '\\Log.txt'),status='replace')
	close(100)
!
!	Read input data
	call LecDat()
!
!	Read sampler data
	call LecSampler()
!
!	Read topo
	write(6,'("------------------------------------")')
	write(6,'("Read dem.asc...",/)')
	call LecTopo()  
!
!	Backup initial DEM values
	topo = topoIni    
!
!	Fill sinks process
	write(6,'("Geometry preprocessing...",/)')
	write(6,'("Fill sinks...")')
	call FillSinksCalc()
!
!	Compute slopes
	write(6,'("Compute slopes...")')
	call SlopeCalc()
!
!	Compute cum area
	write(6,'("Compute flow accumulation...",/)')
	call CumFlowCalc(topo, cumflow, Dinf)
!
!	Zones read
	write(6,'("------------------------------------")')
	write(6,'("Read soils.asc...",/)')
	call LecZones()
!
!	Read soil data
	write(6,'("Read soils data...",/)')
	call LecZonesDat()
!
!	Zones read
	write(6,'("------------------------------------")')
	write(6,'("Read land use lulc.asc...",/)')
	call LecLandUse()
!
!	Read soil data
	write(6,'("Read land use data...",/)')
	call LecLandUsesDat()
!
!	Compute Curve Number
	write(6,'("------------------------------------")')
	write(6,'("Compute CN...",/)')    
	call ComputeCN()    
!
!	Compute averaged CN
	write(6,'("Compute weighted CN...",/)')
	call WeightedCumFlowCalc(topo, CNGrid, WeightedCN, Dinf, cumflow)
!
!   Compute soil data Gaussian parameters
	write(6,'("------------------------------------")')
	write(6,'("Compute soils data Gaussian...",/)')
    call GaussianParams()
!
!	Unconditionally instable cells
	write(6,'("------------------------------------")')
	write(6,'("Inconditionally unstable cells calculation...",/)')
	call IncondUnst()
!    
!   Write UU results
    if (iPROB_uncond_unst .EQ. 1) then
        call WriteGrid(UncIns, 'PROB_uncond_unst.asc')
        call Histogram(UncIns, 'PROB_uncond_unst_HIST.csv')
    endif
!
!	Unconditionally stable cells
	write(6,'("------------------------------------")')
	write(6,'("Inconditionally stable cells calculation...",/)')
	call IncondSta()
!    
!   Write US results
    if (iPROB_uncond_stable .EQ. 1) then
        call WriteGrid(UncEst, 'PROB_uncond_stable.asc')
        call Histogram(UncEst, 'PROB_uncond_stable_HIST.csv')
    endif    
!
!   Count no null cells
    numNullCells = COUNT(FS_mu .EQ. nodata)
    numNotNullCells = my*mx - numNullCells
!
!
!   Sampling loop, antecendent rainfall
    do i = 0, iNumAR
!
       	write(6,'("------------------------------------")')
        write(6,'("Antecedent rainfall iteration: ",I8,/)') i
!
        Rainfall_ant = rMinAR +DBLE(i) * (rMaxAR - rMinAR) / DBLE(iNumAR)
!
!       Event rainfall
        do j = 0, iNumER
!            
           	write(6,'("------------------------------------")')
            write(6,'("Antecedent rainfall iteration: ",I8)') i
            write(6,'("Event rainfall iteration: ",I8,/)') j
!
            Rainfall = rMinER +DBLE(j) * (rMaxER - rMinER) / DBLE(iNumER)
!
!	        Compute averaged antedecent rainfall
	        write(6,'("Compute weighted antecedent rainfall...",/)')
	        call WeightedCumFlowCalc(topo, Rainfall_ant, WeightedRainfall_ant, Dinf, cumflow)
!
!	        Compute averaged event rainfall
	        write(6,'("Compute weighted rainfall...",/)')
	        call WeightedCumFlowCalc(topo, Rainfall, WeightedRainfall, Dinf, cumflow)
!
!           Compute infiltration rainfall
	        write(6,'("------------------------------------")')
            write(6,'("Compute infiltrated rainfall...",/)')
            call Hydrology()
!
!           Antecedent rainfall failure probability
	        write(6,'("------------------------------------")')
	        write(6,'("Antecent rainfall condition...",/)')
	        call InitialSaturation()
!
!           Postevent failure probability
	        write(6,'("------------------------------------")')
	        write(6,'("After event stability computation...",/)')
	        call FinalSaturation()
!
!           Store the sampling results
!
!           Convert nodata cells to 0
            PFGrid = -(PFGrid .NE. nodata) * PFGrid
!
!           Compute cells over threshold
            samplingRes(i+1,j+1) = COUNT( PFGrid >= rThres ) / numNotNullCells
!                       
!    
        enddo
    enddo
!
!
!   Write sampling results
    open(unit=100, file=(trim(fname_res) // '\' // 'samplig_res.txt'), ACTION="write", STATUS="replace")
!
!   Sampling loop, antecendent rainfall
    do i = 0, iNumAR
        write(100, '(1000F14.7)')( real(samplingRes(i+1,j+1)) ,j=0,iNumER)
    enddo
!
    close(100)
!
!
!   Free space required for runoff calculations
	DEALLOCATE(PFGrid)
    DEALLOCATE(h_z)
    DEALLOCATE(h_wt)
    DEALLOCATE(Infiltration)
    DEALLOCATE(WeightedRainfall_ant)
	DEALLOCATE(Rainfall_ant)    
	DEALLOCATE(topoIni)
!
!
!	Liberamos memoria
	DEALLOCATE(topo)
	DEALLOCATE(cumflow)
	DEALLOCATE(zones)
    DEALLOCATE(lulc)
	DEALLOCATE(slopeGrid)
    DEALLOCATE(Soils)
    DEALLOCATE(LandUses)
	DEALLOCATE(GaussKs)
	DEALLOCATE(GaussC)
    DEALLOCATE(GaussCr)
	DEALLOCATE(GaussTanPhi)
	DEALLOCATE(Gaussh)
	DEALLOCATE(GaussPor)
	DEALLOCATE(GaussDens)
	DEALLOCATE(Rainfall)
	DEALLOCATE(CNGrid)    
    DEALLOCATE(WeightedRainfall)    
    DEALLOCATE(WeightedCN)
    DEALLOCATE(samplingRes)
!
!
!
end program fslam_sampling

