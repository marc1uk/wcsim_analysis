		if(cursor->ChildCount()){ // if particles associated to parents
			for(size_t iCh = 0; iCh<ds->GetMC()->GetMCParticleCount(); iCh++){ // Loop on all particles in that event
				cursor->GoChild(iCh); // go to particle
				node = cursor->Here(); // "attach" node to this particle
				if( std::find(interest_volumes_mu.begin(), interest_volumes_mu.end(), node->GetVolume()) != interest_volumes_mu.end() ) { // only interaction in volumes of interest
					if(iCh == 0) { Ninteractions_tot++;} // only fill this once per event
					if (node->GetParticleName() == "mu-" || node->GetParticleName() == "mu+") { // if muon
						is_mu_tag = true;
						Nmuons_tot++;
						for(size_t jCh = 0; jCh<cursor->StepCount(); jCh++){ //loop on each step
							node = cursor->GoStep(jCh); // go to step
							if ( std::find(interest_volumes_mu.begin(), interest_volumes_mu.end(), node->GetVolume()) != interest_volumes_mu.end() ) { // is node is in the volume you want
								vMuTrack.push_back(node->GetEndpoint()); // fills the vectors of node position (front() and back() are first and last node in volumes (track) )
								vMuTrack_Edep.push_back(node->GetKE()); // record current KE of muon at each step of the track
							}
						}
						if (vMuTrack.size() != 0){
							Nmuons_track++;
							muTrack_start = vMuTrack.front();
							muTrack_end = vMuTrack.back();
							hTrackLength_mu->Fill((muTrack_end - muTrack_start).Mag());
							hEdep_muTrack->Fill(vMuTrack_Edep.front() - vMuTrack_Edep.back()); // deposited Edep of track if KEfinal-KEinitial
							//				 cout << "Muon track length: " << (muTrack_end - muTrack_start).Mag() << endl;
							if ((muTrack_end - muTrack_start).Mag() > cut_mu_track){ //muon track length cut
								Nmuons_cut++; is_cut_mu_track = true;
							}
						} else {
							//					 cout << "No muon track in this volume\n";
						}
					}
					
				}
				if (node->GetParticleName() == "neutron") {
					pparticles_trackID.push_back(node->GetTrackID());				
				}
				cursor->GoParent(); // "detach" node
			}
		}
		
		
		/////////////////////////////////////////////////////////
		//********** Secondary neutron search loop ************//
		/////////////////////////////////////////////////////////
		//		 for(size_t iCh = 0; iCh<ds->GetMC()->GetMCTrackCount(); ++iCh){ // Loop on all tracks in that event
		while(node != 0){
			node = cursor->FindNextTrack(); // go to the track
			if(node == NULL){break;} // break the loop if no more non-electron tracks
			
			//			 cout << ds->GetMC()->GetMCTrack(iCh)->GetParticleName() << " " << ds->GetMC()->GetMCTrack(iCh)->GetLastMCTrackStep()->GetProcess()	<< endl;
			//			 if (ds->GetMC()->GetMCTrack(iCh)->GetPDGCode() == 2112){ //look for neutrons, all neutrons
			//				 Nneutrons_tot++;
			//				 if (ds->GetMC()->GetMCTrack(iCh)->GetLastMCTrackStep()->GetProcess() == "nCapture" /*&& 
			//					 std::find(interest_volumes.begin(), interest_volumes.end(), ds->GetMC()->GetMCTrack(iCh)->GetLastMCTrackStep()->GetVolume()) != interest_volumes.end() &&
			//					 is_mu_tag*/) {
			//					 cout << ds->GetMC()->GetMCTrack(iCh+1)->GetPDGCode() << " " << ds->GetMC()->GetMCTrack(iCh+1)->GetLastMCTrackStep()->GetProcess() << endl;
			//					 Nneutrons_cut++;
			//					 NCaptures_perevt++; // nb of ncaptures after muon
			//					 nCapture_pos = ds->GetMC()->GetMCTrack(iCh)->GetLastMCTrackStep()->GetEndpoint();
			//				 //						 cout << "Neutron capture at " << nCapture_pos.x() << " " <<	nCapture_pos.y() << " " << nCapture_pos.z() << endl;
			//					 distance_nCap_muTrack = ((nCapture_pos - muTrack_start).Cross(nCapture_pos - muTrack_end)).Mag()/(muTrack_end - muTrack_start).Mag();
			//				 //						 cout << "Distance nCapture to MuonTrack = " << distance_nCap_muTrack << endl;
			//					 hDist_nCap_muTrack->Fill(distance_nCap_muTrack);
			//					 }
			//					 if (ds->GetMC()->GetMCTrack(iCh)->GetLastMCTrackStep()->GetProcess() == "neutronInelastic") {
			// //						 N_inelastic++;
			//					 }
			//			 }
			if (node->GetParticleName() == "neutron"){ // loop on all neutron tracks
				Nneutrons_track_tot++;
			}
			
			//						cout << node->GetParticleName() << " " << node->GetPDGCode() << " " << node->GetVolume() << " " << node->GetProcess() << " " << node->GetKE() << " " << is_mu_tag << endl; 
			
			if (node->GetProcess() == "nCapture" ) { // capture				
				Nneutrons_cap_tot++;
				if( std::find(interest_volumes_neu.begin(), interest_volumes_neu.end(), node->GetVolume()) != interest_volumes_neu.end() ) { // in the good volumes
					Nneutrons_cap_vol++;
					if (is_mu_tag) { // with a tagged muon
						Nneutrons_cap_mu++;
						if (is_cut_mu_track) { // with a tagged muon having a track longer than the threshold cut
						
						if (parenttrackID != cursor->Parent()->GetTrackID() ) {
							if (is_nGd) { // if next neutron as parent, fill info related to this neutron
								parenttrackID = cursor->Parent()->GetTrackID();
								hEdep_muTrack_nCap->Fill(Edep_capture);
								if (Edep_capture > cut_cap_edep) {
									is_cut_cap_edep = true;
								}
								if(is_cut_cap_edep && is_cut_mu_cap_DR && is_cut_mu_cap_DT) {
									Nneutrons_cap_allcut++;
								}
								if (is_cut_cap_edep && is_cut_mu_cap_DT) {
									Nneutrons_cap_DT++;
								}
								if (is_cut_cap_edep) {
									Nneutrons_cap_Ecut++;
								}
//								 cout << "fill edep a: " << Edep_capture << endl;
								Edep_capture = 0; is_nGd = false; is_cut_cap_edep = false; is_cut_mu_cap_DR = false; is_cut_mu_cap_DT = false; 
							} else {
								Edep_capture = 0; is_nGd = false; is_cut_cap_edep = false; is_cut_mu_cap_DR = false; is_cut_mu_cap_DT = false; 
							}
						}
						
						if (TPMERegexp("1000[0-9][0-9][0-9][0-9][0-9][0-9]").Match(Form("%d",node->GetPDGCode()))) { // capture on an atom (no gammas from nCapture))
							is_nGd = false; is_nH = false;
							//							cout << "atom: " << node->GetPDGCode() << " " << node->GetTrackID() << endl;
							//						 cout << "atom parent: " << cursor->Parent()->GetPDGCode() << " " << cursor->Parent()->GetTrackID() << endl;
							if (TPMERegexp("100064[0-9][0-9][0-9][0-9]").Match(Form("%d",node->GetPDGCode()))) { // PDGCode (Int) casted as a TString, look for Gd
								is_nGd = true;
								Nneutrons_cap_gd++;
								Ncaptures_perevt++; // nb of ncaptures after muon (per evt)
								if(std::find(pparticles_trackID.begin(), pparticles_trackID.end(), cursor->Parent()->GetTrackID()) != pparticles_trackID.end()) { // if the parent neutron is a primary particle
									Npneutrons_cap_gd++;
									Npcaptures_perevt++; // nb of primary ncaptures after muon 
								}
								nCapture_pos = node->GetEndpoint();
								distance_nCap_muTrack = ((nCapture_pos-muTrack_start).Cross(nCapture_pos-muTrack_end)).Mag()/(muTrack_end - muTrack_start).Mag();
								if (distance_nCap_muTrack < cut_mu_cap_DR){
									is_cut_mu_cap_DR = true;
								}
								hDist_nCap_muTrack->Fill(distance_nCap_muTrack);
								if (node->GetGlobalTime() < cut_mu_cap_DT){
									is_cut_mu_cap_DT = true;
								}
								hTime_nCap_muTrack->Fill(node->GetGlobalTime());
								
							} else if (TPMERegexp("100001[0-9][0-9][0-9][0-9]").Match(Form("%d",node->GetPDGCode()))) { //	look for H
								is_nH = true;
								// for now do nothing if n-H
							} else {
								// nothing either
							}
						} else if (node->GetParticleName() == "gamma"){
							//							cout << "gamma: " << node->GetPDGCode() << " " << node->GetTrackID() << " " << node->GetKE() << endl;
							//						 cout << "gamma parent: " << cursor->Parent()->GetPDGCode() << " " << cursor->Parent()->GetTrackID() << endl;
							for( Int_t iStep = 1; iStep < cursor->StepCount(); iStep++ ){
								if(std::find(interest_volumes_neuEdep.begin(), interest_volumes_neuEdep.end(), cursor->Step(iStep)->GetVolume()) != interest_volumes_neuEdep.end()){
									Edep_capture += cursor->Step(iStep-1)->GetKE() - cursor->Step(iStep)->GetKE();
									//								 cout << "stepE: " << cursor->Step(iStep-1)->GetKE() - cursor->Step(iStep)->GetKE() << endl;
								}
							}
							
						}
						parenttrackID = cursor->Parent()->GetTrackID();
					}				
				}
				}
			}
		}
