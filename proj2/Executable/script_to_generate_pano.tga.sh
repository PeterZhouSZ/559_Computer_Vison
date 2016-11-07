# convert wall_000.jpg wall_000.tga
# convert wall_001.jpg wall_001.tga
# convert wall_002.jpg wall_002.tga
# convert wall_003.jpg wall_003.tga
# convert wall_004.jpg wall_004.tga
# convert wall_005.jpg wall_005.tga
# convert wall_006.jpg wall_006.tga
# convert wall_007.jpg wall_007.tga
# convert wall_008.jpg wall_008.tga
# convert wall_009.jpg wall_009.tga
# convert wall_010.jpg wall_010.tga
# convert wall_011.jpg wall_011.tga
# convert wall_012.jpg wall_012.tga
# convert wall_013.jpg wall_013.tga
# convert wall_014.jpg wall_014.tga
# convert wall_015.jpg wall_015.tga


# ./Panorama sphrWarp wall_000.tga wall_000_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_001.tga wall_001_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_002.tga wall_002_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_003.tga wall_003_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_004.tga wall_004_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_005.tga wall_005_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_006.tga wall_006_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_007.tga wall_007_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_008.tga wall_008_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_009.tga wall_009_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_010.tga wall_010_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_011.tga wall_011_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_012.tga wall_012_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_013.tga wall_013_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_014.tga wall_014_warp.tga 571 0.15 0.0
# ./Panorama sphrWarp wall_015.tga wall_015_warp.tga 571 0.15 0.0

# convert wall_000_warp.tga wall_000_warp.pgm
# convert wall_001_warp.tga wall_001_warp.pgm
# convert wall_002_warp.tga wall_002_warp.pgm
# convert wall_003_warp.tga wall_003_warp.pgm
# convert wall_004_warp.tga wall_004_warp.pgm
# convert wall_005_warp.tga wall_005_warp.pgm
# convert wall_006_warp.tga wall_006_warp.pgm
# convert wall_007_warp.tga wall_007_warp.pgm
# convert wall_008_warp.tga wall_008_warp.pgm
# convert wall_009_warp.tga wall_009_warp.pgm
# convert wall_010_warp.tga wall_010_warp.pgm
# convert wall_011_warp.tga wall_011_warp.pgm
# convert wall_012_warp.tga wall_012_warp.pgm
# convert wall_013_warp.tga wall_013_warp.pgm
# convert wall_014_warp.tga wall_014_warp.pgm
# convert wall_015_warp.tga wall_015_warp.pgm

# ./sift < campus_000_warp.pgm > campus_000_warp.key
# ./sift < campus_001_warp.pgm > campus_001_warp.key
# ./sift < campus_002_warp.pgm > campus_002_warp.key
# ./sift < campus_003_warp.pgm > campus_003_warp.key
# ./sift < campus_004_warp.pgm > campus_004_warp.key
# ./sift < campus_005_warp.pgm > campus_005_warp.key
# ./sift < campus_006_warp.pgm > campus_006_warp.key
# ./sift < campus_007_warp.pgm > campus_007_warp.key
# ./sift < campus_008_warp.pgm > campus_008_warp.key
# ./sift < campus_009_warp.pgm > campus_009_warp.key
# ./sift < campus_010_warp.pgm > campus_010_warp.key
# ./sift < campus_011_warp.pgm > campus_011_warp.key
# ./sift < campus_012_warp.pgm > campus_012_warp.key
# ./sift < campus_013_warp.pgm > campus_013_warp.key
# ./sift < campus_014_warp.pgm > campus_014_warp.key
# ./sift < campus_015_warp.pgm > campus_015_warp.key
# ./sift < campus_016_warp.pgm > campus_016_warp.key
# ./sift < campus_017_warp.pgm > campus_017_warp.key

./Features matchSIFTFeatures wall_000_warp.key wall_001_warp.key 0.8 000-001.match 2
./Features matchSIFTFeatures wall_001_warp.key wall_002_warp.key 0.8 001-002.match 2
./Features matchSIFTFeatures wall_002_warp.key wall_003_warp.key 0.8 002-003.match 2
./Features matchSIFTFeatures wall_003_warp.key wall_004_warp.key 0.8 003-004.match 2
./Features matchSIFTFeatures wall_004_warp.key wall_005_warp.key 0.8 004-005.match 2
./Features matchSIFTFeatures wall_005_warp.key wall_006_warp.key 0.8 005-006.match 2
./Features matchSIFTFeatures wall_006_warp.key wall_007_warp.key 0.8 006-007.match 2
./Features matchSIFTFeatures wall_007_warp.key wall_008_warp.key 0.8 007-008.match 2
./Features matchSIFTFeatures wall_008_warp.key wall_009_warp.key 0.8 008-009.match 2
./Features matchSIFTFeatures wall_009_warp.key wall_010_warp.key 0.8 009-010.match 2
./Features matchSIFTFeatures wall_010_warp.key wall_011_warp.key 0.8 010-011.match 2
./Features matchSIFTFeatures wall_011_warp.key wall_012_warp.key 0.8 011-012.match 2
./Features matchSIFTFeatures wall_012_warp.key wall_013_warp.key 0.8 012-013.match 2
./Features matchSIFTFeatures wall_013_warp.key wall_014_warp.key 0.8 013-014.match 2
./Features matchSIFTFeatures wall_014_warp.key wall_015_warp.key 0.8 014-015.match 2
./Features matchSIFTFeatures wall_015_warp.key wall_000_warp.key 0.8 015-000.match 2


./Panorama alignPair wall_000_warp.key wall_001_warp.key 000-001.match 200 1 sift
./Panorama alignPair wall_001_warp.key wall_002_warp.key 001-002.match 200 1 sift
./Panorama alignPair wall_002_warp.key wall_003_warp.key 002-003.match 200 1 sift
./Panorama alignPair wall_003_warp.key wall_004_warp.key 003-004.match 200 1 sift
./Panorama alignPair wall_004_warp.key wall_005_warp.key 004-005.match 200 1 sift
./Panorama alignPair wall_005_warp.key wall_006_warp.key 005-006.match 200 1 sift
./Panorama alignPair wall_006_warp.key wall_007_warp.key 006-007.match 200 1 sift
./Panorama alignPair wall_007_warp.key wall_008_warp.key 007-008.match 200 1 sift
./Panorama alignPair wall_008_warp.key wall_009_warp.key 008-009.match 200 1 sift
./Panorama alignPair wall_009_warp.key wall_010_warp.key 009-010.match 200 1 sift
./Panorama alignPair wall_010_warp.key wall_011_warp.key 010-011.match 200 1 sift
./Panorama alignPair wall_011_warp.key wall_012_warp.key 011-012.match 200 1 sift
./Panorama alignPair wall_012_warp.key wall_013_warp.key 012-013.match 200 1 sift
./Panorama alignPair wall_013_warp.key wall_014_warp.key 013-014.match 200 1 sift
./Panorama alignPair wall_014_warp.key wall_015_warp.key 014-015.match 200 1 sift
./Panorama alignPair wall_015_warp.key wall_000_warp.key 015-000.match 200 1 sift



# Create pairlist.txt from the above results.

# ./Panorama blendPairs pairlist.txt pano.tga 100
