const API_URL = 'http://localhost:3000/api';

let images = [];
let currentIndex = 0;
let cropper = null;
let rotation = 0;
let currentDirectory = '';

// DOM elements
const imageDisplay = document.getElementById('imageDisplay');
const photoDirInput = document.getElementById('photoDir');
const setDirBtn = document.getElementById('setDirBtn');
const loadBtn = document.getElementById('loadBtn');
const prevBtn = document.getElementById('prevBtn');
const nextBtn = document.getElementById('nextBtn');
const rotateLeftBtn = document.getElementById('rotateLeftBtn');
const rotateRightBtn = document.getElementById('rotateRightBtn');
const cropBtn = document.getElementById('cropBtn');
const resetCropBtn = document.getElementById('resetCropBtn');
const saveBtn = document.getElementById('saveBtn');
const resetBtn = document.getElementById('resetBtn');
const currentImageSpan = document.getElementById('currentImage');
const totalImagesSpan = document.getElementById('totalImages');
const imageListContainer = document.getElementById('imageListContainer');
const statusDiv = document.getElementById('status');

// Initialize
document.addEventListener('DOMContentLoaded', () => {
  loadImages();
  setupEventListeners();
});

function setupEventListeners() {
  setDirBtn.addEventListener('click', setDirectory);
  loadBtn.addEventListener('click', loadImages);
  prevBtn.addEventListener('click', previousImage);
  nextBtn.addEventListener('click', nextImage);
  rotateLeftBtn.addEventListener('click', () => rotateImage(-90));
  rotateRightBtn.addEventListener('click', () => rotateImage(90));
  cropBtn.addEventListener('click', toggleCrop);
  resetCropBtn.addEventListener('click', resetCrop);
  saveBtn.addEventListener('click', saveImage);
  resetBtn.addEventListener('click', resetImage);

  // Keyboard shortcuts
  document.addEventListener('keydown', (e) => {
    if (e.target.tagName === 'INPUT') return;

    switch(e.key) {
      case 'ArrowLeft':
        previousImage();
        break;
      case 'ArrowRight':
        nextImage();
        break;
      case 's':
        if (e.ctrlKey || e.metaKey) {
          e.preventDefault();
          saveImage();
        }
        break;
    }
  });
}

async function setDirectory() {
  const directory = photoDirInput.value.trim();
  if (!directory) {
    showStatus('Please enter a directory path', 'error');
    return;
  }

  try {
    const response = await fetch(`${API_URL}/set-directory`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ directory })
    });

    const data = await response.json();

    if (data.success) {
      currentDirectory = data.directory;
      showStatus(`Directory set to: ${data.directory}`, 'success');
      loadImages();
    } else {
      showStatus(data.error || 'Failed to set directory', 'error');
    }
  } catch (error) {
    showStatus('Error setting directory: ' + error.message, 'error');
  }
}

async function loadImages() {
  try {
    showStatus('Loading images...', 'info');
    const response = await fetch(`${API_URL}/images`);
    const data = await response.json();

    images = data.images;
    currentDirectory = data.directory;
    photoDirInput.value = currentDirectory;

    if (images.length === 0) {
      showStatus('No JPEG images found in directory', 'error');
      totalImagesSpan.textContent = '0';
      currentImageSpan.textContent = '0';
      return;
    }

    totalImagesSpan.textContent = images.length;
    currentIndex = 0;
    renderImageList();
    loadCurrentImage();
    showStatus(`Loaded ${images.length} images from ${currentDirectory}`, 'success');
  } catch (error) {
    showStatus('Error loading images: ' + error.message, 'error');
  }
}

function renderImageList() {
  imageListContainer.innerHTML = '';
  images.forEach((image, index) => {
    const item = document.createElement('div');
    item.className = 'image-list-item';
    item.textContent = image;
    item.addEventListener('click', () => {
      currentIndex = index;
      loadCurrentImage();
      renderImageList();
    });

    if (index === currentIndex) {
      item.classList.add('active');
    }

    imageListContainer.appendChild(item);
  });
}

function loadCurrentImage() {
  if (images.length === 0) return;

  destroyCropper();
  rotation = 0;

  const currentImage = images[currentIndex];
  imageDisplay.src = `${API_URL}/image/${encodeURIComponent(currentImage)}`;

  currentImageSpan.textContent = currentIndex + 1;
  updateButtonStates();
  renderImageList();

  showStatus(`Viewing: ${currentImage}`, 'info');
}

function previousImage() {
  if (currentIndex > 0) {
    currentIndex--;
    loadCurrentImage();
  }
}

function nextImage() {
  if (currentIndex < images.length - 1) {
    currentIndex++;
    loadCurrentImage();
  }
}

function rotateImage(degrees) {
  rotation += degrees;
  rotation = rotation % 360;

  if (cropper) {
    cropper.rotate(degrees);
  } else {
    imageDisplay.style.transform = `rotate(${rotation}deg)`;
  }

  showStatus(`Rotated ${degrees > 0 ? 'right' : 'left'}`, 'info');
}

function toggleCrop() {
  if (cropper) {
    destroyCropper();
    cropBtn.textContent = 'Enable Crop';
    cropBtn.classList.remove('btn-warning');
    cropBtn.classList.add('btn-primary');
    showStatus('Crop mode disabled', 'info');
  } else {
    imageDisplay.style.transform = '';
    cropper = new Cropper(imageDisplay, {
      viewMode: 1,
      dragMode: 'move',
      aspectRatio: NaN,
      autoCropArea: 1,
      restore: false,
      guides: true,
      center: true,
      highlight: false,
      cropBoxMovable: true,
      cropBoxResizable: true,
      toggleDragModeOnDblclick: false,
    });
    cropBtn.textContent = 'Disable Crop';
    cropBtn.classList.remove('btn-primary');
    cropBtn.classList.add('btn-warning');
    showStatus('Crop mode enabled - drag to crop', 'info');
  }
}

function resetCrop() {
  if (cropper) {
    cropper.reset();
    showStatus('Crop reset', 'info');
  }
}

function destroyCropper() {
  if (cropper) {
    cropper.destroy();
    cropper = null;
  }
}

function resetImage() {
  rotation = 0;
  destroyCropper();
  loadCurrentImage();
  showStatus('Image reset', 'info');
}

async function saveImage() {
  if (images.length === 0) return;

  try {
    showStatus('Saving image...', 'info');

    let imageData;

    if (cropper) {
      // Get cropped image
      const canvas = cropper.getCroppedCanvas();
      imageData = canvas.toDataURL('image/jpeg', 0.95);
    } else if (rotation !== 0) {
      // Create canvas for rotation
      const canvas = document.createElement('canvas');
      const ctx = canvas.getContext('2d');
      const img = new Image();

      await new Promise((resolve, reject) => {
        img.onload = resolve;
        img.onerror = reject;
        img.src = imageDisplay.src;
      });

      // Handle rotation
      const radians = rotation * Math.PI / 180;
      const sin = Math.abs(Math.sin(radians));
      const cos = Math.abs(Math.cos(radians));

      canvas.width = img.width * cos + img.height * sin;
      canvas.height = img.width * sin + img.height * cos;

      ctx.translate(canvas.width / 2, canvas.height / 2);
      ctx.rotate(radians);
      ctx.drawImage(img, -img.width / 2, -img.height / 2);

      imageData = canvas.toDataURL('image/jpeg', 0.95);
    } else {
      showStatus('No changes to save', 'info');
      return;
    }

    const currentImage = images[currentIndex];

    const response = await fetch(`${API_URL}/save`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        filename: currentImage,
        imageData: imageData,
        rotation: rotation
      })
    });

    const data = await response.json();

    if (data.success) {
      showStatus(`Saved: ${currentImage}`, 'success');
      // Reload the image to show saved version
      setTimeout(() => loadCurrentImage(), 500);
    } else {
      showStatus('Error saving: ' + data.error, 'error');
    }
  } catch (error) {
    showStatus('Error saving image: ' + error.message, 'error');
  }
}

function updateButtonStates() {
  prevBtn.disabled = currentIndex === 0;
  nextBtn.disabled = currentIndex >= images.length - 1;
}

function showStatus(message, type = 'info') {
  statusDiv.textContent = message;
  statusDiv.className = `status ${type}`;
}
