function H = BarPlot(X)
% =======================================================================
% Creates a bar graph with positive data values stacked on the positive
% quadrant and negative data values stacked on the negative quadrant
% =======================================================================
% H = BarPlot(X)
% -----------------------------------------------------------------------
% INPUT
%   - X: data to plot [nobs x nvars]
% -----------------------------------------------------------------------
% OUTPUT
%   - H: handle to graph
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

H1 = bar((X).*(X>0),'stacked'); 
% for ii=1:length(H1)
%     H1(ii).EdgeColor = H1(ii).FaceColor;
% end
hold on;
H2 = bar((X).*(X<0),'stacked');
for ii=1:length(H2)
    H2(ii).FaceColor = H1(ii).FaceColor;
    H2(ii).EdgeColor = H1(ii).EdgeColor;
end
