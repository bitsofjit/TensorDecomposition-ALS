 function plotHydrogen(eVecs, X, Y, Z, n)
 % INPUT
 % eVecs = the eigenvectors to plot (each eigenvector is a column of eVecs)
 % X,Y,Z = the meshgrids from hydrogen3dFD()
 % n = the number of grid points in each spatial direction
 
    N = size(eVecs,2);
 
    % Create a figure and axes
    f = figure('Visible','off');
    ax = axes('Units','pixels');
    
    % Create slider
    eFunc = uicontrol('Style', 'slider',...
        'Min',1,'Max',N,'Value',1,'SliderStep',[1/(N-1) 1/(N-1)],...
        'Position', [20 355 100 20],...
        'Callback', @changeVec);        

    % Create slider
    contourSlider = uicontrol('Style', 'slider',...
        'Min',.002,'Max',.05,'Value',.02,'SliderStep',[1/30 1/15],...
        'Position', [400 355 120 20],...
        'Callback', @changeVec); 
					
    % text for eigenvector level
    txt = uicontrol('Style','text',...
        'Position',[20 380 120 20],...
        'String','Eigenvector');
    
    % text for contour level
    txt2 = uicontrol('Style','text',...
        'Position',[400 380 120 20],...
        'String','Contour level');
    
    % Make figure visble after adding all components
    changeVec;
    view(3)
    f.Visible = 'on';
    % This code uses dot notation to set properties. 
    % Dot notation runs in R2014b and later.
    % For R2014a and earlier: set(f,'Visible','on');


     function changeVec(source, callbackdata)
         val = eFunc.Value;
         cla
         title(horzcat('Eigenvector: ',num2str(val)));
         isosurface(X,Y,Z,reshape(abs(eVecs(:,val)),n,n,n),contourSlider.Value);
         
         axis equal
         xlim([min(min(min(X))) max(max(max(X)))]);
         ylim([min(min(min(Y))) max(max(max(Y)))]);
         zlim([min(min(min(Z))) max(max(max(Z)))]);
%          Lplot = 20;
%          xlim([-Lplot Lplot]);
%          ylim([-Lplot Lplot]);
%          zlim([-Lplot Lplot]);
         camlight('headlight');
     end
end