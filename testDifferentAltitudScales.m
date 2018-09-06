[x,y,z] = sphere(16); % generate sample data

% scatter3(x(:),y(:),z(:),'bo') % plot sample data
scatter3(x(:),y(:),100*z(:),'bo') % plot sample data

hAx2 = axes('Position',get(gca,'Position')); %create duplicate axes at same position

scatter3(hAx2,x(:),y(:),z(:),'r*') %plot scaled data on duplicate axes
set(hAx2,'Visible','off'); %make duplicate axes invisible

%SO IN SUMMARY : THE AXES REMAINING ARE THE 1ST SO WE HAVE TO PLOT THE
%BIGGEST DATASS IN FIRSTTTTTTTTTTT