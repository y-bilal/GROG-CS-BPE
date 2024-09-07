function [kx,ky,dcf]=make_radial_traj(proj,N,os);
% proj = 50; os= 1;
% N = 100;
% For ESMRMB Parallel Imaging Workshop 2010 in Wuerzburg, Germany
% Written by Nicole Seiberlich on 06/23/2010
% Please contact nes30@case.edu with comments, questions, suggestions.

step=1/os;
ang=[0:180/proj:180-180/proj];
li=-N/2:step:(N/2)-step;
ky=li'*sin(ang.*pi/180);
kx=li'*cos(ang.*pi/180);
[psize,rsize]=size(kx);
dcf=abs(li)'*ones(1,rsize);

% % plot of trajectory spoke by spoke
% 
% figure();hold on
% for(kkk = 1:size(k,2))
%     
% stem(kx(:,kkk),ky(:,kkk),'filled', 'MarkerSize' , 1, 'LineStyle','none');
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
%  pause(1);
%     F(kkk) = getframe(gcf) ;
%       drawnow
% end
% hold off
% % create the video writer with 1 fps
%   writerObj = VideoWriter('radial_animation.avi');
%   writerObj.FrameRate = 1;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

