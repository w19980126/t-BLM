function begin = TriggerTime(data,t,Fs)
[~, begin.Pike] = max(diff(data(:,1)));
[~, begin.Cs] = max(diff(data(:,2)));
if begin.Pike >= begin.Cs
    begin.Frame = 0;
else
    begin.Frame = (begin.Cs - begin.Pike)/10000*Fs;
end
end