function img = ex_7_sevilmis(img, filter)

    for i = 2 : size(img,1) - 1
       
        for j = 2 : size(img,2) - 1
           
           pixel = img(i, j);
           new_pixel = round(pixel);
           img(i, j) = new_pixel;
           err = pixel - new_pixel;
           
           if (size(filter,1) == 3) 
               img(i, j + 1) = img(i, j + 1) + err * (7/16);
               img(i + 1, j) = img(i+1, j) + err * (5/16);
               img(i + 1, j - 1) = img(i + 1, j - 1) + err * (3/16);
               img(i + 1, j + 1) = img(i + 1, j + 1) + err * (1/16);
           elseif (size(filter,1) == 5)&&(i < size(img, 1)-1 && j < size(img, 2)-1)&& j>2
               img(i, j + 1) = img(i, j + 1) + err * (8/42);
               img(i + 1, j) = img(i + 1, j) + err * (8/42);
               img(i, j + 2) = img(i, j + 2) + err * (4/42);
               img(i + 2, j) = img(i + 2, j) + err * (4/42);
               img(i + 1, j - 1) = img(i + 1, j - 1) + err * (4/42);
               img(i + 1, j + 1) = img(i + 1, j + 1) + err * (4/42);
               img(i + 1, j - 2) = img(i + 1, j - 2) + err * (2/42);
               img(i + 1, j + 2) = img(i + 1, j + 2) + err * (2/42);
               img(i + 2, j - 1) = img(i + 2, j - 1) + err * (2/42);
               img(i + 2, j + 1) = img(i + 2, j + 1) + err * (2/42);
               img(i + 2, j - 2) = img(i + 2, j - 2) + err * (1/42);
               img(i + 2, j + 2) = img(i + 2, j + 2) + err * (1/42);
           end
           
        end
    end
end