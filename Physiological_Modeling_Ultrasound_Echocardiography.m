% BMENE6003 Homework 01
% Michael Yang
% February 1, 2021


% Load .AVI files into matlab
Apical = 'Apical4chVolunteer.avi';
Mitral = 'SA_MitralVolunteer.avi';
Papilary = 'SA_PapilaryVolunteer.avi';
Apex = 'SA_ApexVolunteer.avi';

Apic_vid = VideoReader(Apical);
Mit_vid = VideoReader(Mitral);
Pap_vid = VideoReader(Papilary);
Apex_vid = VideoReader(Apex);

% Based on the scale on the side of the image, where 15 cm is max, and 5 cm
% = 120 pixels (measured with ruler), the coversion factor:
conversion = 1/24;
% Taken as every 24 pixels is 1 cm

Samples = [1 2 3 4];

for i = Samples
    if i == 1 %Meaning Apical
        % Apical End Diastole long-axis
        play11 = implay(Apical);
        set(play11.Parent, 'Name', 'Apical 4 Chambers - Record Frame Number for End Diastole');
        waitfor(play11);
        prompt1 = {'Enter Frame Number:'};
        dlgtitle = 'Choose Frame for Apical End Diastole';
        answerApicalED = inputdlg(prompt1, dlgtitle);
        % Extract the individual number of the frame for processing
        answerApicalED_num = answerApicalED{1};
        ans11 = str2double(answerApicalED_num);
        this_frame11 = read(Apic_vid, ans11);
        
        % Use imtool to select the measurement
        tool11 = imtool(this_frame11);
        waitfor(tool11);
        prompt11 = {'Enter ED long-axis distance:'};
        dlgtitle = 'Input Measured Long-Axis Distance';
        input11 = inputdlg(prompt11, dlgtitle);
        input11_num = input11{1};
        ED_ans = str2double(input11_num);
        ED_cm = ED_ans * conversion;

        % Apical End-Systole Long-Axis
        play12 = implay(Apical);
        set(play12.Parent, 'Name', 'Apical 4 Chambers - Record End Systole Frame');
        waitfor(play12);
        prompt2 = {'Enter Frame Number:'};
        dlgtitle2 = 'Choose Frame for Apical End Systole';
        answerApicalES = inputdlg(prompt2, dlgtitle2);
        answerApicalES_num = answerApicalES{1};
        ans12 = str2double(answerApicalES_num);
        this_frame12 = read(Apic_vid, ans12);
        
        tool12 = imtool(this_frame12);
        waitfor(tool12);
        prompt12 = {'Enter ES long-axis distance:'};
        dlgtitle = 'Input Measured long-Axis Distance';
        input12 = inputdlg(prompt12, dlgtitle);
        input12_num = input12{1};
        ES_ans = str2double(input12_num);
        ES_cm = ES_ans * conversion;
        
        fig = uifigure;
        selection = uiconfirm(fig, 'Moving Onto Mitral Video', 'Inter-frame notice');
        
    elseif i == 2 %Meaning Mitral
        % Mitral End Diastole short-axis
        play21 = implay(Mitral);
        set(play21.Parent, 'Name', 'Mitral - Record End-Diastole frame');
        waitfor(play21);
        prompt3 = {'Enter Frame Number:'};
        dlgtitle3 = 'Choose Frame for Mitral End Diastole';
        answerMitralED = inputdlg(prompt3, dlgtitle3);
        answerMitralED_num = answerMitralED{1};
        ans21 = str2double(answerMitralED_num);
        this_frame21 = read(Mit_vid, ans21);
        
        % Use imtool to select the measurement of previous
        tool21 = imtool(this_frame21);
        waitfor(tool21);
        prompt21 = {'Enter ED short-axis distance:'};
        dlgtitle = 'Input Measured Short-Axis Distance';
        input21 = inputdlg(prompt21, dlgtitle);
        input21_num = input21{1};
        SA_ans_MitED = str2double(input21_num);
        SA_cm_MitED = SA_ans_MitED * conversion;
        
        % Mitral End Systole short-axis
        play22 = implay(Mitral);
        set(play22.Parent, 'Name', 'Mitral - Record End-Systole frame');
        waitfor(play22);
        prompt4 = {'Enter Frame Number:'};
        dlgtitle4 = 'Choose Frame for Mitral End Systole';
        answerMitralES = inputdlg(prompt4, dlgtitle4);
        answerMitralES_num = answerMitralES{1};
        ans22 = str2double(answerMitralES_num);
        this_frame22 = read(Mit_vid, ans22);
        
        % Use imtool to select the measurement of previous
        tool22 = imtool(this_frame22);
        waitfor(tool22);
        prompt22 = {'Enter ES short-axis distance:'};
        dlgtitle = 'Input Measured Short-Axis Distance';
        input22 = inputdlg(prompt22, dlgtitle);
        input22_num = input22{1};
        SA_ans_MitES = str2double(input22_num);
        SA_cm_MitES = SA_ans_MitES * conversion;
        
        fig = uifigure;
        selection = uiconfirm(fig, 'Moving Onto Papilary Video', 'Inter-frame notice');
        
    elseif i == 3 %Meaning Papilary
        % Papilary End Diastole short-axis
        play31 = implay(Papilary);
        set(play31.Parent, 'Name', 'Papilary - Record End-Diastole frame');
        waitfor(play31);
        prompt5 = {'Enter Frame Number:'};
        dlgtitle5 = 'Choose Frame for Papilary End Diastole';
        answerPapED = inputdlg(prompt5, dlgtitle5);
        answerPapED_num = answerPapED{1};
        ans31 = str2double(answerPapED_num);
        this_frame31 = read(Pap_vid, ans31);
        
        % Use imtool to select the measurement of previous
        tool31 = imtool(this_frame31);
        waitfor(tool31);
        prompt31 = {'Enter ED short-axis distance:'};
        dlgtitle = 'Input Measured Short-Axis Distance';
        input31 = inputdlg(prompt31, dlgtitle);
        input31_num = input31{1};
        SA_ans_PapED = str2double(input31_num);
        SA_cm_PapED = SA_ans_PapED * conversion;
        
        % Papilary End Systole short-axis
        play32 = implay(Papilary);
        set(play32.Parent, 'Name', 'Papilary - Record End-Systole frame');
        waitfor(play32);
        prompt6 = {'Enter Frame Number:'};
        dlgtitle6 = 'Choose Frame for Papilary End Systole';
        answerPapES = inputdlg(prompt6, dlgtitle6);
        answerPapES_num = answerPapES{1};
        ans32 = str2double(answerPapES_num);
        this_frame32 = read(Pap_vid, ans32);
        
        % Use imtool to select the measurement of previous
        tool32 = imtool(this_frame32);
        waitfor(tool32);
        prompt32 = {'Enter ES short-axis distance:'};
        dlgtitle = 'Input Measured Short-Axis Distance';
        input32 = inputdlg(prompt32, dlgtitle);
        input32_num = input32{1};
        SA_ans_PapES = str2double(input32_num);
        SA_cm_PapES = SA_ans_PapES * conversion;
        
        fig = uifigure;
        selection = uiconfirm(fig, 'Moving Onto Apex Video', 'Inter-frame notice');
        
    elseif i == 4
        % Apex End Diastole short-axis
        play41 = implay(Apex);
        set(play41.Parent, 'Name', 'Apex - Record End-Diastole frame');
        waitfor(play41);
        prompt7 = {'Enter Frame Number:'};
        dlgtitle7 = 'Choose Frame for Apex End Diastole';
        answerApexED = inputdlg(prompt7, dlgtitle7);
        answerApexED_num = answerApexED{1};
        ans41 = str2double(answerPapED_num);
        this_frame41 = read(Apex_vid, ans41);
        
        % Use imtool to select the measurement of previous
        tool41 = imtool(this_frame41);
        waitfor(tool41);
        prompt41 = {'Enter ED short-axis distance:'};
        dlgtitle = 'Input Measured Short-Axis Distance';
        input41 = inputdlg(prompt41, dlgtitle);
        input41_num = input41{1};
        SA_ans_ApexED = str2double(input41_num);
        SA_cm_ApexED = SA_ans_ApexED * conversion;

        % Apex End Systole short-axis
        play42 = implay(Apex);
        set(play42.Parent, 'Name', 'Apex - Record End-Systole frame');
        waitfor(play42);
        prompt8 = {'Enter Frame Number:'};
        dlgtitle8 = 'Choose Frame for Apex End Systole';
        answerApexES = inputdlg(prompt8, dlgtitle8);
        answerApexES_num = answerApexES{1};
        ans42 = str2double(answerPapES_num);
        this_frame42 = read(Apex_vid, ans42);
        
        % Use imtool to select the points of previous
        tool42 = imtool(this_frame42);
        waitfor(tool42);
        prompt42 = {'Enter ES short-axis distance:'};
        dlgtitle = 'Input Measured Short-Axis Distance';
        input42 = inputdlg(prompt42, dlgtitle);
        input42_num = input42{1};
        SA_ans_ApexES = str2double(input42_num);
        SA_cm_ApexES = SA_ans_ApexES * conversion;
        
        fig = uifigure;
        selection = uiconfirm(fig, 'Moving To Display Values', 'Notice');
       
    end
end

% Calculations

% End-Diastole Frames, Simpson's Rule
% By treating the obtained values as diameter
A1_ED = pi * (SA_cm_MitED/2)^2;
A2_ED = pi * (SA_cm_PapED/2)^2;
A3_ED = pi * (SA_cm_ApexED/2)^2;
h_ED = (ED_cm / 3);
EDV = ((A1_ED + A2_ED)*h_ED) + ((A3_ED * h_ED)/2) + ((pi * (h_ED^3)) / 6);

% End-Systole Frames, Simpson's Rule
A1_ES = pi * (SA_cm_MitES/2)^2;
A2_ES = pi * (SA_cm_PapES/2)^2;
A3_ES = pi * (SA_cm_ApexES/2)^2;
h_ES = (ES_cm / 3);
ESV = ((A1_ES + A2_ES)*h_ES)+ ((A3_ES * h_ES)/2) + ((pi * (h_ES^3)) / 6);

% Stroke Volume
SV = EDV - ESV;

% Ejection Fraction
EF = SV / EDV;

% Cardiac Output, in mL/min
HeartRate = 60;
CO = SV * HeartRate;

% Creating message box that will display calculations to user
formatspec1 = 'The End Diastolic Ventricular Cavity Volume is: %d';
msg1 = sprintf(formatspec1, EDV);
formatspec2 = 'The End Systolic Ventricular Cavity Volume is: %d';
msg2 = sprintf(formatspec2, ESV);
formatspec3 = 'The Stroke Volume is: %d';
msg3 = sprintf(formatspec3, SV);
formatspec4 = 'The Ejection Fraction is: %d';
msg4 = sprintf(formatspec4, EF);
formatspec5 = 'The Cardiac Output in mL/min given HR=60, is: %d';
msg5 = sprintf(formatspec5, CO);
f = msgbox({msg1, msg2, msg3, msg4, msg5}, 'Calculated Values');
