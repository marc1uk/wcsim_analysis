/* vim:set noexpandtab tabstop=4 wrap */
void cMRDSubEvent::Print(){                     // a generic print

	// Raw Info not printed:
//	std::vector<Int_t> digi_numphots;           // number of true photons for each digit
//	std::vector<Double_t> digi_phot_ts;         // true hit times of photons in a digit
//	std::vector<Int_t> digi_phot_parents;       // wcsim track IDs of parents that provided photons for a digit
cout<<"NEW SUBEVENT"<<endl
	<<"wcsimfile: "<<wcsimfile<<endl
	<<"run: "<<run_id<<", event: "<<event_id<<", trigger: "<<trigger<<", subevent: "<<mrdsubevent_id<<endl
	<<"num reconstructed tracks: "<<tracksthissubevent.size()<<endl
	<<"num true tracks: "<<truetracks.size()<<endl
	<<"num digits: "<<digi_ids.size()<<endl;
	int digitsintracks=0;
	for(auto&& thetrack : tracksthissubevent) digitsintracks+=thetrack.digi_ids.size();
cout<<"number of digits in tracks: "<<digitsintracks<<" ("<<digi_ids.size()-digitsintracks<<" remain)"<<endl
	<<"num pmts hit: "<<pmts_hit.size()<<endl
	<<"digit times: ";
	for(auto atime : digi_ts) cout<<atime<<", ";
cout<<endl<<"digit charges: ";
	for(auto acharge : digi_qs) cout<<acharge<<", ";
cout<<endl<<"layers hit: ";
	for(auto alayer : layers_hit) cout<<alayer<<", ";
cout<<endl<<"energy deposited: ";
	for(auto aqdeposit : eDepsInLayers) cout<<aqdeposit<<", ";
	cout<<endl;
}
