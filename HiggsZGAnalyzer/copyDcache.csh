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

if ($#argv == 4) then
  set dir  = $1
  set username = $2
  set dirout  = $3
  set usernameOut = $4
endif

ls /eos/uscms/store/user/$username/$dir/*.root  > DataFiles.txt

if (! -d /pnfs/cms/WAX/11/store/user/$usernameOut/$dirout) then
     mkdir -p /pnfs/cms/WAX/11/store/user/$usernameOut/$dirout
endif

set count = 1
set nFiles = `cat DataFiles.txt | wc -l`

while($count <= $nFiles)
  set filepath = `cat DataFiles.txt | head -n $count | tail -1`
  set filename = `cat DataFiles.txt | head -n $count | tail -1 | cut -d "/" -f 9`
  echo "Copying /eos/uscms/store/user/$username/$dir/$filename to /pnfs/cms/WAX/11/store/user/$usernameOut/$dirout/"
  dccp $filepath /pnfs/cms/WAX/11/store/user/$usernameOut/$dirout/.  
  @ count ++
end 

rm DataFiles.txt
