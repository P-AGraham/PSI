copylines = 0; nblocks = 200; nlines = 200; if copylines == nil or nblocks == nil or nlines == nil then io.stderr:write("failed to parse arguments"); os.exit(2); end; yvaries = 'x'; yvaries = (yvaries == 'y'); infile = io.open("MAIN_contourtmp0.dat", "r"); outfile = io.open("MAIN_contourtmp0.table", "w"); mesh = PrepcMesh.new(yvaries, nblocks, nlines, copylines, infile, outfile); corners = false; mesh:autocontour(5, -24.9987335000000000, 25.000000000, corners); io.close(infile); mesh:printcontours(); io.close(outfile);