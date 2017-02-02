<template>
  <div id="app">
    <!-- TopNav-->
    <top-nav class="layout-nav"></top-nav>

    <!-- Page -->
    <div class="layout-page">
      <router-view></router-view>
    </div>

  </div>
</template>

<script>
  import TopNav from './components/TopNav'
  import { bus } from './bus.js'

  export default {
    name: 'app',
    components: {
      'top-nav': TopNav,
    }
  }
</script>

<style>
  #app {
    font-family: 'Avenir', Helvetica, Arial, sans-serif;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
    color: #2c3e50;
  }

  body {
    margin: 0px;
  }

  .layout-nav {
    top: 0px;
    width: 100%;
    position: fixed;
    z-index: 10;
  }

  .layout-page {
    padding-top: 60px;
  }

</style>