MODULE mod_calendar
   
   USE mod_precdef
   USE mod_time
   USE mod_param
   
   INTEGER, DIMENSION(10000,12)       :: daysInMonth   
   
   CONTAINS
     
     SUBROUTINE init_calendar
     ! --------------------------------------------------
     ! 
     ! Purpose
     ! -------
     ! Initialize calendar in TRACMASS
     ! Should be called in readfield if ints == intstart
     ! and also each time step
     ! 
     ! Method
     ! ------
     ! Populate the daysInMonth array 
     ! Set currDay,currHour etc to startDay,startHour 
     !
     ! -------------------------------------------------- 
     
        INTEGER              :: jyear,iyear                
        
        IF (.not. noleap) THEN
           DO jyear=1,10000
              iyear = startYear + jyear - 1
              daysInMonth(jyear,:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
              IF ( MOD(iyear,4) == 0 ) THEN
                 daysInMonth(jyear,2) = 29
                 IF ( MOD(iyear,100) == 0 .AND. MOD(iyear,400) /= 0 ) THEN
                    daysInMonth(jyear,2) = 28
                 END IF
              END IF           
           END DO
        ELSE 
           DO jyear=1,10000
              daysInMonth(jyear,:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
           END DO
        END IF
        
        currSec = startSec
        currMin = startMin
        currHour = startHour
        currDay = startDay
        currMon = startMon
        currYear = startYear
        
     RETURN              
     END SUBROUTINE init_calendar
     
     SUBROUTINE update_calendar
     ! ---------------------------------------------------
     ! 
     ! Purpose
     ! -------
     !
     ! Update date and time of TRACMASS 
     !
     ! Method
     ! ------
     !  
     ! Add some minutes to the clock and update dates accordingly
     !
     ! ---------------------------------------------------
     
     INTEGER                 :: iyear
     
     iyear = currYear - startYear + 1 
     
     ! Find number of minutes to add
     IF (ngcm_unit == 1) THEN
        currStep = ngcm
     ELSE IF (ngcm_unit == 2) THEN
        currStep = ngcm * 60
     ELSE IF (ngcm_unit == 3) THEN
        currStep = ngcm * 60 * 60
     ELSE IF (ngcm_unit == 4) THEN
        currStep = ngcm * 24 * 60 * 60
     ELSE IF (ngcm_unit == 5) THEN
        currStep = daysInMonth(iyear,currMon) * 24 * 60 * 60
     ELSE IF (ngcm_unit == 6) THEN
        currStep = SUM(daysInMonth(iyear,:)) * 24 * 60 * 60
     ELSE
        PRINT*," The ngcm_unit ",ngcm_unit," is not recognised "
        PRINT*," Valid units are 1 (second), 2 (minute), 3 (hour) "
        PRINT*," 4 (day), 5 (month), 6 (year) "
        PRINT*," You can also code new units in calendar.f95"
        PRINT*," For now, I stop "
        STOP
     END IF
     
     ! Now update the time and date
     currSec  = currSec + currStep * nff
     
     ! If currSec > 60 we update the minutes, 
     ! and hours etc until we have currSec < 60 again
     DO WHILE (currSec >= 60)
        currSec = currSec - 60
        currMin = currMin + 1
        IF (currMin >= 60) THEN
           currMin = currMin - 60 
           currHour = currHour + 1
           IF (currHour >= 24) THEN
              currHour = currHour - 24
              currDay  = currDay + 1
              IF (currDay > daysInMonth(iyear,currMon)) THEN
                 currDay = currDay - daysInMonth(iyear,currMon)
                 currMon = currMon + 1
              END IF
              IF (currMon > 12) THEN
                 currMon = currMon - 12
                 currYear = currYear + 1
              END IF
           END IF
        END IF
     END DO          
     
     IF (loopYears) THEN
        IF (currYear > loopEndYear) THEN
           currYear = loopStartYear
           loopIndex = loopIndex + 1
        END IF
     END IF
     
     ! If currSec < 0 (backward trajs) we update the minutes,
     ! and hours etc until we have currSec >= 0 again 
     DO WHILE (currSec < 0)
        currSec = currSec + 60
        currMin = currMin - 1
        IF (currMin < 0) THEN
           currMin = currMin + 60
           currHour = currHour - 1
           IF (currHour < 0) THEN
              currHour = currHour + 24
              currDay  = currDay - 1
              IF (currDay < 0) THEN
                 currMon = currMon - 1
                 currDay = currDay + daysInMonth(iyear,currMon)
              END IF
              IF (currMon < 0) THEN
                 currMon = currMon + 12
                 currYear = currYear - 1
              END IF
           END IF
        END IF
     END DO
     
     RETURN     
     END SUBROUTINE update_calendar
     
END MODULE mod_calendar
