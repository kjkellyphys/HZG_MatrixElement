
 if ( removeMcPhoJ ){
  if( (Abs(phoGenMomPID[iPho])==11  || Abs(phoGenMomPID[iPho])==13  || Abs(phoGenMomPID[iPho])==15)
   && (Abs(phoGenGMomPID[iPho])==11 || Abs(phoGenGMomPID[iPho])==13 || Abs(phoGenGMomPID[iPho])==15))   continue;
  if( (Abs(phoGenMomPID[iPho])==11 || Abs(phoGenMomPID[iPho])==13 || Abs(phoGenMomPID[iPho])==15)
   && (Abs(phoGenGMomPID[iPho])==23 || Abs(phoGenGMomPID[iPho])==24  ))   continue;
  if ( phoGenMomPID[iPho]==22 && Abs(phoGenGMomPID[iPho]==11)) continue;
  if ( phoGenMomPID[iPho]==21 || Abs(phoGenMomPID[iPho])<7)   continue;
  } // MC bracket
