function MorphologyBack = flipmorph(Morphology)

    MorphologyBack = Morphology;

    for i=1:length(Morphology)
        switch Morphology{i}
            case 'Flat'
                MorphologyBack{i} = 'Flat';
            case 'RandomUpright'
                MorphologyBack{i} = 'RandomUpright';
            case 'Inverted'
                MorphologyBack{i} = 'RegularUpright';
            case 'RegularUpright'
                MorphologyBack{i} = 'Inverted';
            case 'RegularUpright24'                     % special stuff below
                MorphologyBack{i} = 'Inverted24';
            case 'RegularUpright36'
                MorphologyBack{i} = 'Inverted36';
        end
    end


end

