#
# Check the CRCs for files shared with photo proper, as listed in crc.dat
#
proc check_photo_crc {{crc_file "crc.dat"} {print 0}} {
   if {![file readable $crc_file]} {
      error "Cannot open $crc_file"
   }
   
   set fd [open $crc_file r]
   while {[gets $fd line] != -1} {
      regexp {^[ 	]*([^ 	]+)([ 	]+([^ 	]*))?} $line foo file foo crc
      
      if {$file == "#"} {
	 if $print {
	    puts $line
	 }
	 
	 continue;
      }

      regexp {([^/]+)$} $file foo lfile
      set crc [crcCalcFromFile $lfile]
      set ncrc [crcCalcFromFile $file]

      if {$crc != $ncrc} {
	 lappend errs [list $file $crc $ncrc]
      }

      if $print {
	 echo [format "%-20s\t0x%x" $file $ncrc]
      }
   }

   if [info exists errs] {
      echo "\nTEST-ERR: Errors in checksums:"

      foreach err $errs {
	 echo [format "TEST-ERR: %-20s\t0x%x (should be: 0x%x)" \
		   [lindex $err 0] [lindex $err 1] [lindex $err 2]]
      }
      return 1
   }

   return 0
}

if [check_photo_crc] {
   flush stdout
   error "Found some changed files"
}
