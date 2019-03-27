function ECLDATA = Eclipse(ECLDATA, nwp)
    
    if(isempty(ECLDATA))
        ECLDATA.Time = 0;   %Total production prediction time
        ECLDATA.CompLoc = [30, 30, 1, 1];  %well compeletion location i, j ,k
        ECLDATA.WHLoc = 30;
    end
    
    %now we want to change the well location and compleletion in the data file
    for i = 1:size(ECLDATA.CompLoc(:,1))
        newLoc(i,:) = ECLDATA.CompLoc(i,:);
    end
    
    %Now we should open the sched file
    file = fopen('sched.dat', 'r');
    C = textscan(file, '%s', 'Delimiter', '');
    fclose(file);
    C = C{:};
    
    %Finding welspecs and compdat line
    La = ~cellfun(@isempty, strfind(C, 'WELSPECS'));    %Creates a vector of 0 and 1 that 1 represents the presence of WELSPECS
    
    %Now we should determine the line number of WELSPECS and COMPADT
    wLN = find(La);
    
    %changing the well location and completeion
    for i = 1:size(ECLDATA.CompLoc(:,1))
            
        %And the next line of COMPDAT is the well compeletion data (cLN+i+nwp)
        %Now we have to split the lines containing data
        SplitStr_wellsp = regexp(cell2mat(C(wLN+i+nwp)), '\ ', 'split');
        
        SplitStr_wellsp{3} = num2str(newLoc(i,1));
        SplitStr_wellsp{4} = num2str(newLoc(i,2));
        newWHLoc = strjoin(SplitStr_wellsp);          % Removing Qoutation
        C{wLN+i+nwp} = newWHLoc;
    end
    %Wrting to file
    fopen('sched.dat', 'w');
    fprintf(file, '%s\r\n', C{:});
    fclose(file);    
   
    %now our new file is ready to run
    dos('$multirun.bat');
    
    A = importdata('SPE10.RSM', ' ', 7);
    
    AA = cell2mat(A.textdata(5)); %extracting the multiplier of FGPT (10^3)
    AAA = str2num(AA(47));        %extracting the power of multiplier (3)
    
    %now check if it has got a value
    if(isempty(AAA))
        AAA = 0;
    end
    
    BBB = str2num(AA(60));        %extracting the multiplier for FWPT
    
    %now check if it has got a value
    if(isempty(BBB))
        BBB = 0;
    end
    
    CCC = str2num(AA(99));        %extracting the multiplier for FWPT
    %now check if it has got a value
    if(isempty(CCC))
        CCC = 0;
    end
    
    
    Time = A.data(20:end, 2);           % Get time in years strating from 1
    FOPT = A.data(20:end, 4)*10^AAA;    % Field gas production rate -- FGPR STB
    FWPT = A.data(20:end, 5)*10^BBB;    % Fied Water production rate -- FWPR STB
    FWIT = A.data(20:end, 8)*10^CCC;
    
    i = 1;
    while (i <= numel(Time))
        if(mod(Time(i), 1))
            Time(i) = [];
            FOPT(i) = [];
            FWPT(i) = [];
            FWIT(i) = [];
            i = i - 1;
        end
        i = i + 1;
    end
    
    ECLDATA.Time = Time;
    ECLDATA.Qop = FOPT;
    ECLDATA.Qwp = FWPT;
    ECLDATA.Qwi = FWIT;
end
