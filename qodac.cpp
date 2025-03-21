void __fastcall DeviceQODAC::generateEvents( const double maxTime, FormulaEnvironment* env ) throw( EEventException ) {
	debug("QODAC: generateEvents...");

	uint8_t num_channels = 0;
	uint8_t channel_states = 0;

	std::vector<uint8_t> lin_error_comp(8);
	std::vector<uint32_t> initial_value(8 , 524288);
	std::vector<uint16_t> num_events(8);
	typedef std::vector<uint32_t> vec_u32;
	std::vector<vec_u32> sampling_periods;
	std::vector<vec_u32> num_samples;
	std::vector<uint32_t> samples;

	int idx_channel = 0;
	FormulaVariable* tVar = new FormulaVariable;
	tVar->name = "t";
	for (IOChannelListIterator itC = channels.begin() ; itC != channels.end() ; itC++) {
		QODACChannel* channel = static_cast<QODACChannel*>(*itC);
		if (!channel->enabled || channel->events.empty()) continue;
		int k = channel->number;
		num_channels++;
		channel_states |= (1<<idx_channel);
		if ((channel->hardwareMin == -10) && (channel->hardwareMax == 10))
			lin_error_comp.at(k) = 12;
		else
			lin_error_comp.at(k) = 0;
		channel->prepareFunction();
		initial_value.at(k) = applyFunction(channel, 0, channel->initialValue);
		num_events.at(k) = channel->events.size();
		std::vector<uint32_t> sampling_periods_ch;
		std::vector<uint32_t> num_samples_ch;
		QODACEventFormula* currentEvent = NULL;
		QODACEventFormula* nextEvent = NULL;
		TimingEvent* tNext = NULL;
		double timeCurrentEvent=0, timeNextEvent;
		for (IOEventListIterator itE = channel->events.begin(); itE != channel->events.end(); ++itE) {
			nextEvent = static_cast<QODACEventFormula*>(*itE);
			tNext = nextEvent->getTimingEvent();
			timeNextEvent = roundTime(tNext->getValue());
			if (!tNext->isEnabled()) continue;
			if (currentEvent == NULL) {
				currentEvent = nextEvent;
				timeCurrentEvent = timeNextEvent;
				continue;
			}
			try {
				tVar->setValue(timeCurrentEvent);
				currentEvent->f.environment->addVariable(tVar);
				currentEvent->f.prepareEvaluate();
				if (currentEvent->timeDependent) {
					double dt = 1/(double)currentEvent->eventSamplingRate;
					sampling_periods_ch.push_back(dt*1e6);
					int i = 0;
					for(double t=timeCurrentEvent; t<timeNextEvent; t+=dt) {
						tVar->setValue(t);
						samples.push_back(applyFunction(channel, tVar->value, currentEvent->f.evaluate()));
						i++;
					}
					num_samples_ch.push_back(i);
				} else {
						sampling_periods_ch.push_back((timeNextEvent-timeCurrentEvent)*1e6);
						num_samples_ch.push_back(1);
						samples.push_back(applyFunction(channel, tVar->value, currentEvent->f.evaluate()));
				}
			} catch (Exception& e) {
				throw EEventException( e.Message , currentEvent );
			}
			currentEvent->f.environment->removeVariable(tVar);
			currentEvent = nextEvent;
			timeCurrentEvent = timeNextEvent;
		}
		try {
			timeNextEvent = roundTime(maxTime);
			tVar->setValue(timeCurrentEvent);
			currentEvent->f.environment->addVariable(tVar);
			currentEvent->f.prepareEvaluate();
			if (currentEvent->timeDependent) {
				double dt = 1/(double)currentEvent->eventSamplingRate;
				double samp_per = dt*1e6;
				int samp_per_int = samp_per;
				sampling_periods_ch.push_back(samp_per_int);
				int i = 0;
				for(double t=timeCurrentEvent; t<timeNextEvent; t+=dt) {
					tVar->setValue(t);
					int kek = applyFunction(channel, tVar->value, currentEvent->f.evaluate());
					samples.push_back(applyFunction(channel, tVar->value, currentEvent->f.evaluate()));
					i++;
				}
				num_samples_ch.push_back(i);
			} else {
					sampling_periods_ch.push_back((timeNextEvent-timeCurrentEvent)*1e6);
					num_samples_ch.push_back(1);
					int kek = applyFunction(channel, tVar->value, currentEvent->f.evaluate());
					samples.push_back(applyFunction(channel, tVar->value, currentEvent->f.evaluate()));
			}
		} catch (Exception& e) {
			throw EEventException( e.Message , currentEvent );
		}
		currentEvent->f.environment->removeVariable(tVar);
		channel->finalizeFunction();
		sampling_periods.push_back(sampling_periods_ch);
		num_samples.push_back(num_samples_ch);
		idx_channel++;
	}
	delete tVar;

	uint8_t debug_messages = 0;
	uint8_t command = 1;
	uint8_t play_mode = 1;
	uint8_t triggerExternal = 0;
	uint8_t clockExternal = 0;
	channel_states = 255;
	uint32_t config_flag = debug_messages + (command << 1) + (play_mode << 4) + (triggerExternal << 5) + (clockExternal << 7) + (channel_states << 8);

	uint64_t metadata_size = 8*3;
	for (int i=0; i<num_channels; i++)
		metadata_size += 2*num_events[i];
	buffer = (uint32_t*)malloc(sizeof(uint32_t)*(1+metadata_size+samples.size()));
	buffer[0] = config_flag;
	buffer_idx++;

	for (int i = 0 ; i < 8 ; i++) {
		buffer[buffer_idx] = lin_error_comp[i];
		buffer_idx++;
	}
	for (int i = 0 ; i < 8 ; i++) {
		buffer[buffer_idx] = initial_value[i];
		buffer_idx++;
	}
	for (int i = 0 ; i < 8 ; i++) {
		buffer[buffer_idx] = num_events[i];
		buffer_idx++;
	}

	for (int i=0; i<sampling_periods.size(); i++) {
		for (int j=0; j<sampling_periods[i].size(); j++) {
			buffer[buffer_idx] = sampling_periods[i][j];
			buffer_idx++;
		}
	}
	for (int i = 0 ; i < num_samples.size() ; i++) {
		for (int j = 0 ; j < num_samples[i].size() ; j++) {
			buffer[buffer_idx] = num_samples[i][j];
			buffer_idx++;
		}
	}
	for (unsigned int i = 0 ; i < samples.size() ; i++) {
		buffer[buffer_idx] = samples.at(i);
		buffer_idx++;
	}
	buffer_idx *= 4;

	debug("QODAC: ...generateEvents");
}

