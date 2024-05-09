function output = medListCheck(medListOpt)

% medListDir = 'G:\My Drive\PHCP\SFM\Paper\Translational Psych Submission\Response to Reviews\';
% medListDir = 'E:\UMN Drive\PHCP\SFM\Paper\Translational Psych Submission\Response to Reviews\';
medListDir = '/home/shaw-raid1/data/psychophysics/SFM.git/';


%% Import list of antipsychotic/antidepressants/benzo meds
psychMedTable = readtable(fullfile(medListDir,'AntipsychoticMedicationTypes.xlsx'));
% psychMedTable = table2cell(psychMedTable);
% Make list of strings from table for meds
for iI = 1:height(psychMedTable)
    medList{iI} = psychMedTable{iI,1}{1};
    benzoList(iI) = psychMedTable{iI,2};
end

% Import list of other non-related drugs
otherMedTable = readtable(fullfile(medListDir,'NonPsychMedList.xlsx'));
% otherMedTable = table2cell(otherMedTable);
% Make large list that combines all other med types
medTypes = {'VitaminsAndSupplements', 'Cholesterol', 'Sleep', 'Respitory', 'Allergy', 'Diabetes', ...
    'BloodPressure', 'PainAndMuscleRelaxers', 'BirthControl', 'Other'};
counter = 0;
for iI = 1:length(medTypes)
    for iJ = 1:length(otherMedTable.(medTypes{iI}))
        if ~isempty(otherMedTable.(medTypes{iI}){iJ})
            counter = counter + 1;
            otherMedList{counter} = otherMedTable.(medTypes{iI}){iJ};
        end
    end
end

%% Grab patient med data
patientTable = readtable(fullfile(medListDir,'PHCP-MedicationHistory_DATA_LABELS_2023-11-16_1140.csv'));
% patientTable = table2cell(patientTable);

% Clean up the patient med table
% First, only look at clinical lists, as these are the most accuracte med
% lists.
patientTable = patientTable(strcmp(patientTable{:,2},'Clinical'),:);

% Only look at included subjects
includedPart = zeros([height(patientTable) 1]);
for iI = 1:height(patientTable)
    if sum(strcmp(medListOpt.partName,patientTable{iI,1}))>0
        includedPart(iI) = 1;
    end
end

% Now grab all reported meds for each participant
% Total of 17 possible medication entries
% 1st med is on column 6
% Grab meds by looking every 6 columns for med name and stop if an empty
% cell is found
counter = 0;
for iPart = 1:height(patientTable)
    if includedPart(iPart) == 1
        counter = counter + 1;
        patientMedList.partID{counter} = patientTable{iPart,1}{1};
        patientMedList.partNum(counter) = str2num(patientTable{iPart,1}{1}(2:end));
        patientMedList.medName{counter} = {};
        for iMed = 1:17
            if ~isempty(patientTable{iPart,iMed*6}{1})
                patientMedList.medName{counter}{iMed} = patientTable{iPart,iMed*6}{1};
            else
                break
            end
        end
    end
end

%% Now go through and check each med and whether:
% 1) It's included in our med list at all and,
% 2) it's a benzodiazapine.

% Find meds that aren't included in our psych meds lists or other meds list
counter = 0;
isNotInMedTable = {};
for iPart = 1:length(patientMedList.partNum)
    for iMed = 1:length(patientMedList.medName{iPart})
        if sum(strcmpi(patientMedList.medName{iPart}{iMed}, medList))==0 & ...
                sum(strcmpi(patientMedList.medName{iPart}{iMed}, otherMedList))==0
            % If we haven't already included this med
            if sum(strcmpi(patientMedList.medName{iPart}{iMed}, isNotInMedTable))==0
                counter = counter+1;
                isNotInMedTable{counter} = patientMedList.medName{iPart}{iMed};
            end
        end
    end
end
isNotInMedTable = isNotInMedTable';

% If there are meds that aren't already in our table of known meds then
% stop and check what they are for potential benzos we're missing
if ~isempty(isNotInMedTable)
    % Write this list to excel sheet for processing
    filename = 'excludedMedList.xlsx';
    delete filename
    writecell(isNotInMedTable,filename);
    error('Check excludedMedList for missing meds!')
end

% Now search the patient med lists for meds that are benzodiazapines and
% make an exclusion index for patients taking benzos
patientMedList.isTakingBenzo = zeros([1 length(patientMedList.partNum)]);
for iPart = 1:length(patientMedList.partNum)
    for iMed = 1:length(patientMedList.medName{iPart})
        if sum(strcmpi(patientMedList.medName{iPart}{iMed}, medList(benzoList==1)))==1
            patientMedList.isTakingBenzo(iPart) = patientMedList.isTakingBenzo(iPart)+1;
            patientMedList.isTakingWhichBenzo{iPart}{patientMedList.isTakingBenzo(iPart)} = patientMedList.medName{iPart}{iMed};
        end
    end
    if patientMedList.isTakingBenzo(iPart) == 0
        patientMedList.isTakingWhichBenzo{iPart}= {};
    end
end

% Now make full list (including repeats) by referencing each part in the
% original patient list with the new exclusion criteria
for iPart = 1:length(medListOpt.partNum)
    % Find index for this particpant in exclusion list
    if patientMedList.isTakingBenzo(patientMedList.partNum == medListOpt.partNum(iPart)) > 0
        output(iPart) = 1;
    else
        output(iPart) = 0; 
    end
end

output = output';

end



