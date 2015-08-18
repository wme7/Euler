%demo_textprogressbar
%This a demo for textprogressbar script
textprogressbar('calculating outputs: ',100);
for i=1:100,
    textprogressbar(i,100);
    %pause(0.1);
end
textprogressbar('done');


textprogressbar('saving data:         ',80);
for i=1:0.5:80,
    textprogressbar(i,80);
    %pause(0.05);
end
textprogressbar('terminated');