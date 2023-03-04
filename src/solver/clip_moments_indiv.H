//  Info << "Clip moments" << endl;
  const double SMALLWEIGHT_C = 1.0e-20;
  const double NUCL_NBRC_C   = 2 * 20.0;
  const double NUCL_SURF_C   = std::pow(NUCL_NBRC_C,2.0/3.0);
  // Check for globally small moments
  if ( MomArr[0] < SMALLWEIGHT_C ||
       MomArr[1] < SMALLWEIGHT_C*NUCL_NBRC_C ||
       MomArr[2] < SMALLWEIGHT_C*(std::pow(NUCL_NBRC_C,2)) ||
       MomArr[3] < SMALLWEIGHT_C*(std::pow(NUCL_NBRC_C,3)) ||
       MomArr[4] < SMALLWEIGHT_C*NUCL_SURF_C ||
       MomArr[5] < SMALLWEIGHT_C*NUCL_NBRC_C*NUCL_SURF_C )
  {
     MomArr[1] = max(MomArr[1],SMALLWEIGHT_C*NUCL_NBRC_C);
     MomArr[0] = MomArr[1] / NUCL_NBRC_C;
     MomArr[2] = max(MomArr[2],SMALLWEIGHT_C*(std::pow(NUCL_NBRC_C,2)));
     MomArr[3] = max(MomArr[3],SMALLWEIGHT_C*(std::pow(NUCL_NBRC_C,3)));
     MomArr[4] = NUCL_SURF_C * MomArr[0];
     //MomArr[4] = SMALLWEIGHT_C*NUCL_SURF_C;
     MomArr[5] = max(MomArr[5],SMALLWEIGHT_C*NUCL_NBRC_C*NUCL_SURF_C);
  }

  if (MomArr[1] < NUCL_NBRC_C*MomArr[0] ||
      MomArr[4] < NUCL_SURF_C*MomArr[0])
  {
     MomArr[0] = MomArr[1] / NUCL_NBRC_C;
     MomArr[4] = MomArr[0] * NUCL_SURF_C;
  }

  if (nMom>=6)
  {
     MomArr[2] = max(MomArr[2],MomArr[1]*MomArr[1]/MomArr[0]);
     MomArr[3] = max(MomArr[3],MomArr[1]*MomArr[1]*MomArr[1]/MomArr[0]);
     MomArr[5] = max(MomArr[5],MomArr[1]*MomArr[4]/MomArr[0]);
  }
