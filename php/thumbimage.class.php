<?php
// https://code.tutsplus.com/tutorials/how-to-create-a-thumbnail-image-in-php--cms-36421
class ThumbImage
{
    private $source;
 
    public function __construct($sourceImagePath)
    {
        $this->source = $sourceImagePath;
    }
 
    public function createThumb($destImagePath, $thumbWidth=100)
    {
        $sourceImage = imagecreatefrompng($this->source);
        $orgWidth = imagesx($sourceImage);
        $orgHeight = imagesy($sourceImage);
        $thumbHeight = floor($orgHeight * ($thumbWidth / $orgWidth));
        $destImage = imagecreatetruecolor($thumbWidth, $thumbHeight);
        imagecopyresampled($destImage, $sourceImage, 0, 0, 0, 0, $thumbWidth, $thumbHeight, $orgWidth, $orgHeight);
        imagepng($destImage, $destImagePath);
        imagedestroy($sourceImage);
        imagedestroy($destImage);
    }
}
