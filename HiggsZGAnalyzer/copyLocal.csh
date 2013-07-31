#!/bin/csh

if ($#argv == 2) then
  set dir  = $1
  set username = $2
  set dirout  = $1
endif

if ($#argv == 3) then
  set dir  = $1
  set username = $2
  set dirout  = $3
endif

ls /pnfs/cms/WAX/11/store/user/$username/$dir/*.root  > DataFiles.txt

if (! -d ~/nobackup/$dirout) then
     mkdir -p ~/nobackup/$dirout
endif

set count = 1
set nFiles = `cat DataFiles.txt | wc -l`

while($count <= $nFiles)
  set filepath = `cat DataFiles.txt | head -n $count | tail -1`
  set filename = `cat DataFiles.txt | head -n $count | tail -1 | cut -d "/" -f 11`
  echo "Copying store/user/$username/$dir/$filename to ~/nobackup/$dirout/"
  dccp $filepath ~/nobackup/$dirout/.  
  @ count ++
end 

rm DataFiles.txt
