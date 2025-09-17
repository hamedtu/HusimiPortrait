%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PREPRINT: "Pre-Floquet states facilitating coherent subharmonic   %%%%%
%%%%%%%%%% response of periodically driven many-body systems"   %%%%%%%%%%%%%%%%%%
%%%%%%%%%% Steffen Seligmann, Hamed Koochaki Kelardeh, and Martin Holthaus   %%%%%
%%%%%%%%%%               https://arxiv.org/abs/2504.00578               %%%%%%%%%%
%%%%%%%%%%*******************************************************************%%%%%
%%%%%%%%%%         This script execute datasets and plots Figures      %%%%%%%%%%%
%%%%%%%%%%                  from the above paper                       %%%%%%%%%%%
%%%%%%%%%% For instance  "Fig_choice=1" produces Fig.1, etc. %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all

Fig_choice = 3; %   %%%%%%%%%%  contrl the Figures according to the paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Fig_choice
    case 1
        Mat = readmatrix('Fig1.txt');
    case 2
        Mat = readmatrix('Fig2.txt');
    case 3
        Mat1 = load("TS_eigenvectors_0.92_0.40_1.90_10000_64_real.dat");
        Mat2 = load("TS_eigenvectors_0.92_0.40_1.90_10000_64_imag.dat");
    case 5
        Mat1 = readmatrix('Fig5a.txt');
        Mat2 = readmatrix('Fig5b.txt');
    case 6
        Mat1 = load("TS_eigenvectors_0.92_0.40_1.90_10000_64_real.dat");
        Mat2 = load("TS_eigenvectors_0.92_0.40_1.90_10000_64_imag.dat");
    case 7
        Mat = readmatrix('Fig7.txt');
    case 8
        Mat = readmatrix('Fig8.txt');
    case 9
        Mat1 = readmatrix('Fig9a.txt');
        Mat2 = readmatrix('Fig9b.txt');
    otherwise
        warning('Unexpected scan type. No plot created. ')
end

if Fig_choice==1  % Plotting the Return probability
    tn=Mat(:,1);
    ReturnP=Mat(:,2);
    figure;
    plot(tn,ReturnP,'-k','LineWidth',1)
    set(gca,'ydir','normal','fontsize',24);
    ylabel("| \langle \psi (0)|\psi (t)\rangle | ^2", 'FontSize',28)
    xlabel('t/T','FontSize',28)
    title('Return probability','FontSize',20)
    ylim([0 1])
    xticks([0 3 6 9 12 15])

elseif Fig_choice==2  % Plotting the Poincare map
    ncyc=1000; w= 1.9;
    T=2*pi/w;
    tspan = 0:T:ncyc*T;
    nt = length(tspan);
    figure;
    plot(Mat(nt+1:end,:),Mat(1:nt,:),'.')
    set(gca,'ydir','normal','fontsize',24);
    xlim([-pi pi]);ylim([-1 1])
    xlabel("$\varphi/\pi$",Interpreter="latex",FontSize=28);ylabel("p",'FontSize',28);
    xticks([-pi -pi/2  0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

elseif Fig_choice==3  % Color-coded Husimi projections
    FloquetMatrix = Mat1 + 1i * Mat2;
    N = 10000;
    alpha = 0.92;
    mu = 0.4;
    w = 1.9;

    [eta, I] = Coherence(FloquetMatrix, N);

    % Figure 3
    index = 10001-[1, 110, 194, 276, 768, 972, 1415, 1673];
    Husimiplot(FloquetMatrix, I(index), [alpha mu w], N, 250);
elseif Fig_choice==5
    figure;
    tiledlayout(2,1);
    nexttile
    scatter(Mat1(:,1),Mat1(:,2)./pi,[],Mat1(:,3),".")
    colormap(gca,"jet");
    set(gca,'ydir','normal','fontsize',24);
    ylabel("$\varphi/\pi$",Interpreter="latex",FontSize=28);xlabel("t/T",'FontSize',28);
    axis("tight")
    yticks([-1 -1/2 -1/4  0 1/4 1/2 1])

    nexttile
    scatter(Mat2(:,1),Mat2(:,2)./pi,[],Mat2(:,3),".")
    colormap(gca,"jet");
    set(gca,'ydir','normal','fontsize',24);
    ylabel("$\varphi/\pi$",Interpreter="latex",FontSize=28);xlabel("t/T",'FontSize',28);
    axis("tight")
    yticks([-1 -1/2 -1/4  0 1/4 1/2 1])

elseif Fig_choice==6  % Color-coded Husimi projections
    FloquetMatrix = Mat1 + 1i * Mat2;
    N = 10000;
    alpha = 0.92;
    mu = 0.4;
    w = 1.9;

    [eta, I] = Coherence(FloquetMatrix, N);
    % Figure 6
    index2 = 5705;
    Husimiplot(FloquetMatrix, I(index2), [alpha mu w], N, 250);
elseif Fig_choice==7
    % Plotting the Occupation Probability
    % imagesc('XData',1:HamD,'YData',tn(1:451),'CData',prj_fock(1:451,:) )
    figure;
    imagesc('XData',Mat(:,1),'YData',Mat(1:451,2),'CData',Mat(1:451,3:end) )
    set(gca,'ydir','normal','fontsize',24);
    colormap jet; colorbar;
    axis("tight")
    xlabel('Fock states','FontSize',28,'Color','k' );
    ylabel('t/T ','FontSize',28,'Color','k');
    title('occupation probability','FontSize',20)

elseif Fig_choice==8  % Plotting the Poincare map
    ncyc=999; w= 1.9;
    T=2*pi/w;
    tspan = 0:T:ncyc*T;
    nt = length(tspan);
    figure;
    plot(Mat(nt+1:end,:),Mat(1:nt,:),'.')
    set(gca,'ydir','normal','fontsize',24);
    xlim([-0.1*pi 0.1*pi]);ylim([-0.6 -0.4])
    xlabel("$\varphi/\pi$",Interpreter="latex",FontSize=28);ylabel("p",'FontSize',28);
    xticks([-0.1*pi 0 0.1*pi ])
    xticklabels({'-0.1','0','0.1'})
    yticks([-0.58 -0.5 -0.42 ])
    yticklabels({'-0.58','-0.5','-0.42'})

elseif Fig_choice==9
    figure;
    tiledlayout(2,1);
    nexttile
    % tn=Mat1(:,1);
    % ReturnP=Mat1(:,2);
    plot(Mat1(:,1),Mat1(:,2),'-k','LineWidth',1)
    set(gca,'ydir','normal','fontsize',24);
    ylabel("| \langle \psi (0)|\psi (t)\rangle | ^2", 'FontSize',28)
    xlabel('t/T','FontSize',28)
    title('N=2000','FontSize',20)
    ylim([0 0.6])
    xticks([0 3 6 9 12 15 18 21 24 27 30 33 36 39])

    nexttile
    plot(Mat2(:,1),Mat2(:,2),'-k','LineWidth',1)
    set(gca,'ydir','normal','fontsize',24);
    ylabel("| \langle \psi (0)|\psi (t)\rangle | ^2", 'FontSize',28)
    xlabel('t/T','FontSize',28)
    title('N=5000','FontSize',20)
    ylim([0 0.6])
    xticks([0 3 6 9 12 15 18 21 24 27 30 33 36 39])

end