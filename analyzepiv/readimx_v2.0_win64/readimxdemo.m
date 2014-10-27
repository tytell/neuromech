[lv_pathstr]= fileparts(which('readimx'));
lv_dirlist = dir([lv_pathstr '/TestImages']);
lv_j=1;
clear A;
for lv_i=1:size(lv_dirlist,1),
   if ( ~lv_dirlist(lv_i).isdir)
	   disp(['Load and display ' lv_dirlist(lv_i).name ' ...'] ); 
       A(lv_j) = readimx([lv_pathstr  '/TestImages/' lv_dirlist(lv_i).name]);
       figure; showimx(A(lv_j).Frames{1});title(['Image ' num2str(lv_j) ': ' lv_dirlist(lv_i).name]);
       lv_j = lv_j+1;
   end
end
clear lv_*;