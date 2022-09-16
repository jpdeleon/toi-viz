<div class="container">
<div class="gallery">
  <?php
  require "thumbimage.class.php";

  // Image extensions
  $image_extensions = array("png","jpg","jpeg","gif");

  // Target directory
  $dir = 'images/';
  if (is_dir($dir)){

   if ($dh = opendir($dir)){
    $count = 1;

    // Read files
    while (($file = readdir($dh)) !== false){

     if($file != '' && $file != '.' && $file != '..'){


      // Image path
      $image_path = "images/".$file;
      $objThumbImage = new ThumbImage($image_path);

      // Thumbnail image path
      $thumbnail_path = "images/thumbnail/".$file;

      if (!file_exists($thumbnail_path)) {
          //echo "Creating thumbnail: $thumbnail_path\n";
          $objThumbImage->createThumb($thumbnail_path, 250);
      }
      //else {
      //    echo "Using $thumbnail_path\n";
      //}
      $thumbnail_ext = pathinfo($thumbnail_path, PATHINFO_EXTENSION);
      $image_ext = pathinfo($image_path, PATHINFO_EXTENSION);

      // Check its not folder and it is image file
      if(!is_dir($image_path) &&
         in_array($thumbnail_ext,$image_extensions) &&
         in_array($image_ext,$image_extensions)){
   ?>

       <!-- Image -->
       <a href="<?php echo $image_path; ?>">
        <img src="<?php echo $thumbnail_path; ?>" alt="" title=""/>
       </a>
       <!-- --- -->
       <?php

       // Break
       if( $count%4 == 0){
       ?>
         <div class="clear"></div>
       <?php
       }
       $count++;
      }
     }

    }
    closedir($dh);
   }
  }
 ?>
</div>
</div>
