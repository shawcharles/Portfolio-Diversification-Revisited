function DataOut=CutData(DataIn,pas)
t=1;
count=1;
while t<=length(DataIn)
    DataOut(count)=DataIn(t);
    count=count+1;
    t=t+pas;
end