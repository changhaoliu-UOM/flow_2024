function cmplx = read_tess(filename,cmplx)
% -------------------------------------------------------------------------
% read_tess.m
% -------------------------------------------------------------------------
% open Neper .tess file
% -------------------------------------------------------------------------
fid = fopen(filename);

% -------------------------------------------------------------------------
% get vertex data
% - cmplx(1).num(1).val
% - cmplx(1).bc
% -------------------------------------------------------------------------
while ~feof(fid); if contains(fgetl(fid),cmplx(1).name); break; end; end

% get number of vertices
cmplx(1).num(1).val = str2double(fgetl(fid));

% get raw vertex data and strip
cmplx(1).bc = zeros(cmplx(1).num(1).val,5);
for i=1:cmplx(1).num(1).val; cmplx(1).bc(i,:) = str2num(fgetl(fid)); end
cmplx(1).bc = cmplx(1).bc(:,2:4);

% -------------------------------------------------------------------------
% get edge data
% - cmplx(2).num(1).val
% - cmplx(2).num(2).val
% - cmplx(2).bndop(1).indx
% - cmplx(2).bndop(1).sgn
% -------------------------------------------------------------------------
while ~feof(fid); if contains(fgetl(fid),cmplx(2).name); break; end; end

% get number of edges
cmplx(2).num(2).val = str2double(fgetl(fid));

% get raw edge data and strip
cmplx(2).bndop(1).indx = zeros(cmplx(2).num(2).val,4);
for i=1:cmplx(2).num(2).val
    cmplx(2).bndop(1).indx(i,:) = str2num(fgetl(fid));
end
cmplx(2).bndop(1).indx = cmplx(2).bndop(1).indx(:,2:3);
cmplx(2).bndop(1).sgn = ...
    [ones(cmplx(2).num(2).val,1),-ones(cmplx(2).num(2).val,1)];
cmplx(2).num(1).val=2*ones(cmplx(2).num(2).val,1);

% -------------------------------------------------------------------------
% get face data
% - cmplx(3).num(3).val
% - cmplx(3).num(2).val
% - cmplx(3).bndop(2).indx
% - cmplx(3).bndop(2).sgn
% - cmplx(3).num(1).val      - not needed?
% - cmplx(3).bndop(1).indx   - not needed?
% - cmplx(3).bndop(1).sgn    - not needed?
% - cmplx(3).dir 
% -------------------------------------------------------------------------
while ~feof(fid); if contains(fgetl(fid),cmplx(3).name); break; end; end

% get number of faces
cmplx(3).num(3).val = str2double(fgetl(fid));

% get raw face data and strip
cmplx(3).bndop(1).indx = zeros(cmplx(3).num(3).val,10);
cmplx(3).bndop(2).indx = zeros(cmplx(3).num(3).val,10);
cmplx(3).dir = zeros(cmplx(3).num(3).val,4);
for i=1:cmplx(3).num(3).val
    tmp = str2num(fgetl(fid));
    cmplx(3).bndop(1).indx(i,1:length(tmp)) = tmp;
    tmp = str2num(fgetl(fid));
    cmplx(3).bndop(2).indx(i,1:length(tmp)) = tmp;
    cmplx(3).dir(i,:) = str2num(fgetl(fid));
    junk = str2num(fgetl(fid));
end
cmplx(3).num(1).val=cmplx(3).bndop(1).indx(:,2);
cmplx(3).bndop(1).sgn = ...
    sign(cmplx(3).bndop(1).indx(:,3:max(cmplx(3).num(1).val)+2));
cmplx(3).bndop(1).indx = ...
    abs(cmplx(3).bndop(1).indx(:,3:max(cmplx(3).num(1).val)+2));
cmplx(3).num(2).val=cmplx(3).bndop(2).indx(:,1); 
cmplx(3).bndop(2).sgn = ...
    sign(cmplx(3).bndop(2).indx(:,2:max(cmplx(3).num(2).val)+1));
cmplx(3).bndop(2).indx = ...
    abs(cmplx(3).bndop(2).indx(:,2:max(cmplx(3).num(2).val)+1));
cmplx(3).dir = cmplx(3).dir(:,2:4);

% -------------------------------------------------------------------------
% get poly data
% - cmplx(4).num(4).val
% - cmplx(4).num(3).val
% - cmplx(4).bndop(3).indx
% - cmplx(4).bndop(3).sgn
% -------------------------------------------------------------------------
while ~feof(fid); if contains(fgetl(fid),cmplx(4).name); break; end; end

% get number of polys
cmplx(4).num(4).val = str2double(fgetl(fid));

% get raw poly data and strip
cmplx(4).bndop(3).indx = zeros(cmplx(4).num(4).val,20);
for i=1:cmplx(4).num(4).val
    tmp = str2num(fgetl(fid));
    cmplx(4).bndop(3).indx(i,1:length(tmp)) = tmp;
end
cmplx(4).num(3).val=cmplx(4).bndop(3).indx(:,2); 
cmplx(4).bndop(3).sgn = ...
    sign(cmplx(4).bndop(3).indx(:,3:max(cmplx(4).num(3).val)+2));
cmplx(4).bndop(3).indx = ...
    abs(cmplx(4).bndop(3).indx(:,3:max(cmplx(4).num(3).val)+2));

% get seed locations (Delaunay vertices)
frewind(fid)
while ~feof(fid); if contains(fgetl(fid),'seed'); break; end; end
for i=1:cmplx(4).num(4).val
    tmp = str2num(fgetl(fid)); cmplx(4).cc(i,1:3) = tmp(2:4);
end

% -------------------------------------------------------------------------
% close Neper .tess file
% -------------------------------------------------------------------------
fclose(fid);
