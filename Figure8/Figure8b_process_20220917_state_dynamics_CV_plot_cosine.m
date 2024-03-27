
%{
----------------------------------------------------------------------------
Determine the time window for task relevant neurons
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}
function Figure8b_process_20220917_state_dynamics_CV_plot_cosine

auc_cosine = [82.591
86.816
86.6477
81.3141
89.165
86.4065
86.6828
85.5378
80.7569
55.327
87.0141
79.4032
85.4151
78.101
86.31
78.3189
87.5599
81.2108
78.8685
71.6287
86.5615
88.767
82.5699
88.9869
57.0533
84.1157
77.4173
87.1135
];

fof_cosine = [77.218
40.5717
65.7238
47.7303
65.0712
47.1829
67.9153
47.6635
65.2971
89.6874
21.3405
50.4837
64.7534
77.9371
62.3134
54.4493
81.2793
84.0557
81.8514
86.2847
71.2494
53.7178
41.8168
60.1502
82.807
72.3033
86.4236
87.3553
66.0442
68.914
70.0133
];

mpfc_cosine = [74.1656
62.9424
86.2793
87.9217
86.5414
75.7231
84.1622
66.3885
77.8604
46.5448
74.123
75.3967
29.7308
67.43
85.9532
82.6052
87.6112
89.5988
77.0344
67.4163
64.1827
41.5698
84.8077
87.805
73.2725
87.7152
81.5753
85.1073
82.6496
87.4761
79.6934
78.9448
61.8008
];


ranksum(auc_cosine, fof_cosine)
ranksum(auc_cosine, mpfc_cosine)
ranksum(fof_cosine, mpfc_cosine)

box1 = ones(length(auc_cosine),1);
box2 = ones(length(fof_cosine),1) * 2;
box3 = ones(length(mpfc_cosine),1) * 3;

x1 = (rand(length(auc_cosine),1)-0.5) * 0.2 + 1;
x2 = (rand(length(fof_cosine),1)-0.5) * 0.2 + 2;
x3 = (rand(length(mpfc_cosine),1)-0.5) * 0.2 + 3;

%% Fig7a
cd('G:\upload_code\Figure8\Fig8b');
figure; hold on

boxplot([auc_cosine; fof_cosine; mpfc_cosine], [box1;box2;box3])
plot(x1, auc_cosine, 'k.')
plot(x2, fof_cosine, 'k.')
plot(x3, mpfc_cosine, 'k.')
set(gca,'ylim',[15 90])
sdata.region =[box1;box2;box3];
sdata.cosine = [auc_cosine; fof_cosine; mpfc_cosine];
T = struct2table(sdata);
writetable(T, 'source fig8b.csv');

