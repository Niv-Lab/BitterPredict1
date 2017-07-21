function [bitterPredictOutput] = BitterPredict( propertyFile, varargin)
%BitterPredict is a classifier which predicts weather a given molecule is likely 
%to be bitter or non-bitter. 
%for more information on the method please refer to (link will be available upon
% publication acceptance)
%INPUT: The user should provide xls/xlsx or csv file with calculated properties 
%by Schrodinger physicochemical descriptors and QikProp package. 
% The user can specify additional parameters as detailed in the parameters
% section below.

%The input file must include the following 59 descriptors:
%-------------------------------------------------------------
%{'#stars','#amine','#amidine','#acid','#amide','#rotor','#rtvFG','CNS','dipole','SASA','FOSA','FISA',
%'PISA','volume','donorHB','accptHB','dip^2/V','ACxDN^.5/SA','glob','QPpolrz','QPlogPC16','QPlogPoct',
%'QPlogPw','QPlogPo/w','QPlogS','CIQPlogS','QPlogHERG','QPPCaco','QPlogBB','QPPMDCK','QPlogKp','IP(eV)'
%,'EA(eV)','#metab','QPlogKhsa','HumanOralAbsorption','PercentHumanOralAbsorption','SAfluorine',
%'SAamideO','PSA','RuleOfFive','RuleOfThree','#ringatoms','#in34','#in56','#noncon','#nonHatm','MW',
%'AlogP','RB','HeavyAtomCount','RingCount','Estate','MR','Polar','Num atoms','Num aromatic rings',
%'ChiralCenterCount','Total charge'};
 
%The input file may include other descriptors and/or molecule identifier. % 
%The input file should include header line with the descriptor’s names as appears in the list above.
%Empty/nan values are not allowed in the descriptors listed above.

%----------------------------------------------------------------
%Input File example:  would be available upon publication acceptance via the Git repository 
% **** User can set several parameters*****:
%'modelFile': hold the path to bitterPredict model file (i.e bitterPredict2017.mat file) default:
% = local directory
% 'workSheet': hold excel worksheet name, If the input file is excel user can specify from which 
% worksheet to load data. if not specified the data from the first sheet. 
% 'numberOfIdsColumnsToPrint' –Usually the first column or columns are % molecule identiferes (such as 
% smiles and/or molecules ID)this %parameter determine how many first columns would be printed in the 
% result. For example if numberOfIdsColumnsToPrint'=2 the first two columns would be also printed to 
% the result. Default=1 
%save File: indicates if to save the result file. The result file type will be the same as the input 
% file (csv or excel) default: True. 
%'outputPath' : specify where to save the result file (default = local directory)
%'printStat' : print summary of results (how many  bitter and non-bitter molecules where predicted)
% 'saveResultFile': hold if to save the result file .
%'resultPath': path to result file default = local directory 
%'setName'- if specified the set name will be added to the results file name

%command line examples: 
%xls:
%BitterPredict(‘<Path to input file>\inputFile.xlsx','modelFile', <(‘<Path to model file>\bitterPredict2017.mat','resultPath','<path to folder>');
%csv
%BitterPredict('<Path to input file>\InputData.csv','modelFile', <Path to model file>\bitterPredict2017.mat','resultPath',’<path to result folder> ',’numberOfIdsColumnsToPrint',2);

%% process input
inFile=propertyFile;
propForModelOrdered={'#stars','#amine','#amidine','#acid','#amide','#rotor','#rtvFG','CNS','dipole','SASA','FOSA','FISA','PISA','volume','donorHB','accptHB','dip^2/V','ACxDN^.5/SA','glob','QPpolrz','QPlogPC16','QPlogPoct','QPlogPw','QPlogPo/w','QPlogS','CIQPlogS','QPlogHERG','QPPCaco','QPlogBB','QPPMDCK','QPlogKp','IP(eV)','EA(eV)','#metab','QPlogKhsa','HumanOralAbsorption','PercentHumanOralAbsorption','SAfluorine','SAamideO','PSA','RuleOfFive','RuleOfThree','#ringatoms','#in34','#in56','#noncon','#nonHatm','MW','AlogP','RB','HeavyAtomCount','RingCount','Estate','MR','Polar','Num atoms','Num aromatic rings','ChiralCenterCount','Total charge'};

% available parameters :
bitterPredictParams=getParamsFromUser('modelFile','bitterPredict2017.mat','firstNumericColumn',1 ,'printStat',1,'numberOfIdsColumnsToPrint',0,'saveResultFile',1,'resultPath','.','setName','','workSheet','',varargin{:});
 %load input file
try
    load(bitterPredictParams.modelFile); %model path + file
catch
   display (['ERROR:  failed to load bitterPredict2017.mat.\n  suggested solutions: add bitterPredict2017.mat folder to path, or set the modelFile parmeter to hold full path to bitterPredict2017.mat']);
end

%check if file is excel or csv 
fileNameSplit= strsplit(inFile,'.');
fileExtension=fileNameSplit{end};
fileType='xls';

%% Get input from Excel file 
if(regexp(fileExtension, 'xls*'))
%If worksheet is specified, read from specific worksheet. If not read data from the first worksheet
%in the excel file. 
    try
        if(~isempty(bitterPredictParams.workSheet))
             [num,txt,raw]=xlsread(inFile,bitterPredictParams.workSheet);
        else
            [num,txt,raw]=xlsread(inFile);
        end
    catch
        display ['ERROR read input file TO COMPLETE'];
    end
    %Extract header line from the excelfile
    propHeader=txt(1,:);
    num_cols=size(raw,2);
    num_rows=size(raw,1);
    first_prop=raw{1,1};
    %Get first numeric column in order to know how many first columns of raw cells are not present 
    % in num matrix
    firstNumericColumn=getDataFirstCol(raw);
    bitterPredictParams.firstNumericColumn=max(firstNumericColumn,1);
    %If user did not set a valid number of identifier columns to print all first columns till first 
    % numeric column would be printed. 
     if (( ~isnumeric(bitterPredictParams.numberOfIdsColumnsToPrint)||( bitterPredictParams.numberOfIdsColumnsToPrint)<1 || (bitterPredictParams.numberOfIdsColumnsToPrint>size(raw,2))) )
        bitterPredictParams.numberOfIdsColumnsToPrint=max(bitterPredictParams.firstNumericColumn-1,1);
     end
    %Organize input data to be in the required order for bitterPredict model
    orderedData=[];
    indexsToOrderedProp=zeros(59,1);
    errorWithProp=false;
    for i=1:length(propForModelOrdered)
        currProp=propForModelOrdered{i};
        %check that the current  descriptor exists
        currPropCompareResult=strcmp(lower(propHeader),lower(currProp));
        if(sum(currPropCompareResult)>0)
            indexProp=find(currPropCompareResult,1);
            %copy relevant descriptor values
            orderedData(:,i)=num(:,(indexProp-bitterPredictParams.firstNumericColumn)+1);
            %check it does not contain any nan values
            if(sum(isnan(orderedData(:,i))))
                display(['ERROR: the data include non numeric cells , a non numeric elemtes were found in the : "' currProp '" property please fix the input data']);
                errorWithProp=true;
                return;
            end
            indexsToOrderedProp(i)=indexProp;
        else
            display(['ERROR: property: "', currProp , '" was not found. Please check your input file and make sure it contains the missing property ']);
            errorWithProp=true;
            return;
        end
    end
%% Input is CSV file 
else
    if(strcmp(fileExtension,'csv') || strcmp(fileExtension,'txt'))
        fileType='csv';
        try
        propTable=readtable(inFile);
        catch 
              display ['ERROR reading input file'];
        end 
      
         formatToCSV='';
        num_cols=size(propTable,2);
        num_rows=size(propTable,1)+1;
        first_prop=propTable.Properties.VariableNames{1};
         %If user did not set a valid number of identifier columns to print all first columns till first 
         % numeric column would be printed. 
        if (( ~isnumeric(bitterPredictParams.numberOfIdsColumnsToPrint)||( bitterPredictParams.numberOfIdsColumnsToPrint)<1 || (bitterPredictParams.numberOfIdsColumnsToPrint>size(propTable,2))) )
                    bitterPredictParams.firstNumericColumn=getDataFirstColTable(propTable);
                    bitterPredictParams.numberOfIdsColumnsToPrint=max(bitterPredictParams.firstNumericColumn-1,1);
        end
        orderedData=[];
        errorWithProp=false;
        for i=1:length(propForModelOrdered)
            currProp=propForModelOrdered{i};
            %convert descriptors name to matlab style 
            expressionNumberSign='#';
            replaceWith='x_';
            currProp = regexprep(currProp,expressionNumberSign,replaceWith);
            expressionSpecialChars='[\^\/\.\(\)]';
            replaceWith='_';
            currProp = regexprep(currProp,expressionSpecialChars,replaceWith);
            expressionSpace='\s';
            replaceWith='';
            currProp = regexprep(currProp,expressionSpace,replaceWith);
            %check that all required descriptors exist in the input file
            if(any(strcmp(lower(currProp), lower(propTable.Properties.VariableNames))))
                   propIndex=find(strcmp(lower(currProp), lower(propTable.Properties.VariableNames)));
                   orderedData(:,i)=propTable.(propTable.Properties.VariableNames{propIndex});
            else
                display(['ERROR: property: ' currProp ' was not found. Please check your input file and make sure it contains the missing property ']);
                errorWithProp=true;
                 return;
            end     
            if(sum(isnan(orderedData(:,i))))
                  display(['ERROR: the data include non numeric cells , a non-numeric elemtes were found in the : "' currProp '" property please fix the input data']);
                  errorWithProp=true;
                return;
            end
        end    
    end
end
%% run BitterPrdedict model
if(~errorWithProp)
    %Get only molecules in bitter domain
    [setIndomain,IndexesInBitterDomain,IndexOutOfBitterDomain]=inBitterDomain(orderedData);
    %run BitterPredict AdaBoost and get predictions 
    [predictions,Score]=bitterModel2017.ada200.predict(setIndomain);
    
    %% Prepare output:
    %(in bitter domain default=0, predictions default =nan or 0 , score = -1000)
    bitterPredictOutput=cell(num_rows,(bitterPredictParams.numberOfIdsColumnsToPrint)+3);
    idIsIndexs=false;
    if(bitterPredictParams.numberOfIdsColumnsToPrint==1)
        %check if first prop is part of the data if it is generate an index
        %column to be printed first.
        if strcmp(first_prop, propForModelOrdered)
           bitterPredictOutput(2:end,1)=num2cell((1:num_rows)-1);
           bitterPredictOutput{1,1}='Index';
           idIsIndexs=true;
        end
    end
    if (~idIsIndexs)
       
        for i=1:(bitterPredictParams.numberOfIdsColumnsToPrint) 
            if(fileType=='xls')
                  currNamesCol=raw(:,i);
            else
                currNamesCol=cell(num_rows,1);
                currNamesCol(2:end)=propTable{:,i};
                currNamesCol(1)=propTable.Properties.VariableNames(i);
                formatToCSV=[formatToCSV '%s,'];
            end
             bitterPredictOutput(:,i)=currNamesCol;
        end
    end
    InBitterDomain=zeros(size(orderedData,1),1);
    InBitterDomain(IndexesInBitterDomain)=1;
    out_predictions=nan(size(orderedData,1),1);
    out_predictions(IndexesInBitterDomain)=predictions;
    out_score=nan(size(orderedData,1),1); 
    out_score(IndexesInBitterDomain)=Score(:,2);
    for i=2:size(bitterPredictOutput,1)
        bitterPredictOutput{i,bitterPredictParams.numberOfIdsColumnsToPrint+1}=InBitterDomain(i-1);
        bitterPredictOutput{i,bitterPredictParams.numberOfIdsColumnsToPrint+2}=out_predictions(i-1);
        bitterPredictOutput{i,bitterPredictParams.numberOfIdsColumnsToPrint+3}=out_score(i-1);
    end
    %update headers 
     bitterPredictOutput{1,bitterPredictParams.numberOfIdsColumnsToPrint+1}='InBitterDomain';
     bitterPredictOutput{1,bitterPredictParams.numberOfIdsColumnsToPrint+2}='Prediction';
     bitterPredictOutput{1,bitterPredictParams.numberOfIdsColumnsToPrint+3}='Score';
     
     negativeP=sum(predictions==-1);
     posP=sum(predictions==1);
     if(isfield(bitterPredictParams,'thr'))
        posP=sum(Score(:,2)>=bitterPredictParams.thr);
        negativeP=sum(Score(:,1)<bitterPredictParams.thr);
     end
     if(bitterPredictParams.printStat)
   
        fprintf('Number of compounds ready  prediction: %d\n', num_rows-1);
        fprintf('Number of compounds in "Bitter Domain": %d\n', length(predictions));
        fprintf('Predicted as nonBitter in "Bitter Domain": %d \n',negativeP);
        fprintf('Predicted as Bitter in "Bitter Domain": %d \n',posP);
    
     end
     if(bitterPredictParams.saveResultFile)
         outFile=[bitterPredictParams.resultPath bitterPredictParams.setName 'bitterPredict_'  datestr(now, 'HH_MM_dd-mmm')];
       if(regexp(fileExtension, 'xls*'))
             xlswrite([bitterPredictParams.resultPath bitterPredictParams.setName 'bitterPredict_'  datestr(now, 'HH_MM_dd-mmm')],bitterPredictOutput);
       else
           tmpCell=bitterPredictOutput(2:end,:);
           tableBitterOutput=cell2table(tmpCell);
           tableBitterOutput.Properties.VariableNames = bitterPredictOutput(1,:);
           writetable(tableBitterOutput,[outFile '.txt']);
       end
     end
else
    exit(1);
end

end

function [ outSet,IndexesInBitterDomain,IndexesNotInBitterDomain ] = inBitterDomain( setOrderdProp )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%assuming 59 properties 
MW_index=48;
AlogP_index=49;
MW_thr=700;
AlogP_BottomThr=-3;
AlogPUperThr=7;
IndexesInBitterDomain=(setOrderdProp(:,MW_index)<=MW_thr) & (setOrderdProp(:,AlogP_index)>= AlogP_BottomThr)  & (setOrderdProp(:,AlogP_index)<=AlogPUperThr);
outSet=setOrderdProp(IndexesInBitterDomain,:);
IndexesNotInBitterDomain=~IndexesInBitterDomain;
end

function [predictions_out,Score_out]=AdaBoostPredict(algosStruct)

    [predictions,Score]=algosStruct.bitterPredict.predict(params.outGroupData);
    predictions_out=predictions;
    Score_out=Score;

  
    if(iscell(predictions))
        predictions=cellfun(@str2num,predictions); 
    end

    negativeP=sum(predictions==-1);
    posP=sum(predictions==1);
    
    if(isfield(bitterPredictParams,'thr'))
        posP=sum(Score(:,2)>=bitterPredictParams.thr);
        negativeP=sum(Score(:,1)<bitterPredictParams.thr);
        fprintf('Number of compounds in "Bitter Domain": %d\n', length(predictions));
        fprintf('Predicted as nonBitter above threshold "Bitter Domain": %d \n',negativeP);
        fprintf('Predicted as Bitter in above threshold "Bitter Domain": %d \n',posP);
    end
    
    
end

function [firstDataColumn]=getDataFirstCol(raw)
firstDataRaw=raw(2,:);
numericCol = cellfun(@(x) isnumeric(x) && numel(x)==1, firstDataRaw);
firstDataColumn=find(numericCol,1);

end

function [firstDataColumn]=getDataFirstColTable(propTable)
    firstDataColumn=1;
     for p=1:length(propTable.Properties.VariableNames)
            prop=propTable.Properties.VariableNames{p};
            if (isnumeric(propTable.(prop)))
                firstDataColumn=p;
                break;
            end
     end
end

function [ params ] = getParamsFromUser( varargin )

 numVarargs = length(varargin);
 params=struct();
 j=1;
 while j<numVarargs
    params.(varargin{j})= varargin{j+1};
    j=j+2;
 end
end
 


