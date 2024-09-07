function error=RMSE(org_img,rec_img)
error = sum((org_img(:) - rec_img(:)).^2)/numel(org_img);  % MSE
error = sqrt(error);                                    % RMSE
end
